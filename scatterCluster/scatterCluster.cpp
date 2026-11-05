#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

#include "Model.hpp"
#include "Settings.hpp"
#include "Star.hpp"
#include "Utility.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::fixed;
using std::setprecision;
using std::setw;
using std::string;
using std::vector;

double getExposureTime(const Settings &s, const string &f)
{
    auto it = s.scatterCluster.exposures.find(f);
    if (it != s.scatterCluster.exposures.end())
        return it->second;

    if (s.global.verbose)
        cerr << "No exposure time for filter \"" << f << "\" in scatter config" << endl;

    return 0.0;
}

// Computes S/N based on magnitude and exposure time using coefficients from config
// logS2N = a - b * mag, then scale by sqrt(exptime)
double signalToNoise(const Settings &settings, double mag, double exptime, const string &filter)
{
    auto it = settings.scatterCluster.s2nCoeffs.find(filter);
    if (it == settings.scatterCluster.s2nCoeffs.end())
    {
        cerr << "Warning: No S/N coefficients for filter \"" << filter << "\"" << endl;
        return 0.0;
    }

    double logS2N = it->second.first - it->second.second * mag;
    double s2n = exp10(logS2N);
    s2n *= sqrt(exptime);

    return s2n;
}

bool meetsMagCutoff(const Settings &settings, const StellarSystem &s)
{
    double mag;

    try
    {
        mag = s.obsPhot.at(settings.scatterCluster.relevantFilt);
    }
    catch (std::out_of_range &e)
    {
        cerr << "scatterCluster.relevantFilt out of range" << endl;
        throw;
    }

    if (s.observedStatus == BD)
        return true;
    else if (mag > 99.)
        return false;                       // check if real object.  if not, skip
    else if (mag < settings.scatterCluster.brightLimit)
        return false;                       // typically used to remove red giants
    else if (mag > settings.scatterCluster.faintLimit)
        return false;                       // typically used to remove WDs and lower MS
    else
        return true;
}

bool meetsStageCutoff(const StellarSystem &s, const bool isFS)
{
    if (s.observedStatus == NSBH || s.observedStatus == DNE)
        return false;
    else if (s.observedStatus == WD && s.secondary.mass > 0.0 && !isFS)
        return false;                       // TEMPORARY KLUDGE -- ignore binaries of MS/RG + WDs and WD + WD
    else
        return true;
}


int main(int argc, char *argv[])
{
    Settings settings;

    vector<double> exposureTimes;
    vector<string> filters;
    vector<string> outputFilters;  // Intersection of photometry and config filters
    vector<size_t> filterIndices;  // Indices into original filter list
    vector<StellarSystem> systems;

    settings.loadSettings(argc, argv);

    if (settings.global.seed == std::numeric_limits<uint32_t>::max())
    {
        srand(std::time(0));
        settings.global.seed = rand();

        cout << "Seed: " << settings.global.seed << endl;
    }

    // Determine input photometry file
    string inputFile;
    if (!settings.global.files.phot.empty())
    {
        inputFile = settings.global.files.phot;
    }
    else
    {
        inputFile = settings.global.files.output + ".sim.out";
    }

    {
        // Read photometry
        vector<double> filterPriorMin;
        vector<double> filterPriorMax;

        std::ifstream photin(inputFile);

        if (!photin)
        {
            cerr << "***Error: Photometry file " << inputFile << " was not found.***" << endl;
            cerr << "[Exiting...]" << endl;
            exit(-1);
        }

        auto ret = base::utility::readPhotometry(photin, filterPriorMin, filterPriorMax, settings.global);

        filters = ret.first;

        // Build set of filters that have both exposure times and S/N coefficients in config
        std::unordered_set<string> configFilters;
        for (const auto &kv : settings.scatterCluster.exposures)
        {
            if (kv.second > 0 && settings.scatterCluster.s2nCoeffs.count(kv.first))
            {
                configFilters.insert(kv.first);
            }
        }

        // Compute intersection: filters present in both photometry and config
        for (size_t f = 0; f < filters.size(); ++f)
        {
            if (configFilters.count(filters[f]))
            {
                outputFilters.push_back(filters[f]);
                filterIndices.push_back(f);
                exposureTimes.push_back(getExposureTime(settings, filters[f]));
            }
        }

        if (outputFilters.empty())
        {
            cerr << "No filters in common between photometry and scatter config" << endl;
            cerr << "Photometry filters: ";
            for (const auto &f : filters) cerr << f << " ";
            cerr << endl;
            cerr << "Config filters: ";
            for (const auto &f : configFilters) cerr << f << " ";
            cerr << endl;
            exit(-1);
        }

        if (settings.global.verbose)
        {
            cout << "Using " << outputFilters.size() << " filters (intersection):" << endl;
            for (size_t f = 0; f < outputFilters.size(); ++f)
                cout << "  " << f << " " << outputFilters[f] << " (exposure: " << exposureTimes[f] << "h)" << endl;
        }

        for (auto &s : ret.second)
        {
            bool isFS = std::stoi(s.id) > 1000000;
            if (meetsMagCutoff(settings, s) && meetsStageCutoff(s, isFS))
            {
                bool okay = true;

                for (size_t f = 0; f < outputFilters.size(); ++f)
                {
                    auto filter = outputFilters[f];
                    auto origIdx = filterIndices[f];
                    auto exposureTime = exposureTimes[f];

                    if (exposureTime > 0)
                    {
                        auto s2n = signalToNoise(settings,
                                                 s.obsPhot.at(origIdx),
                                                 exposureTime,
                                                 filter);

                        if (s2n < settings.scatterCluster.limitS2N)
                        {
                            okay = false;
                            break;
                        }
                    }
                }

                if (okay)
                {
                    systems.push_back(s);
                }
            }
        }

        if (systems.empty())
        {
            cerr << "No valid input photometry after considering mag and stage cutoffs" << endl;
            exit(0);
        }
    }

    std::mt19937 gen(uint32_t(settings.global.seed * uint32_t(2654435761))); // Applies Knuth's multiplicative hash for obfuscation (TAOCP Vol. 3)

    for (auto &s : systems)
    {
        for (size_t f = 0; f < outputFilters.size(); ++f)
        {
            auto filter = outputFilters[f];
            auto origIdx = filterIndices[f];
            auto exposureTime = exposureTimes[f];

            if (exposureTime > 0)
            {
                auto s2n = signalToNoise(settings,
                                         s.obsPhot.at(origIdx),
                                         exposureTime,
                                         filter);

                auto sigma = 1 / s2n;

                sigma = std::max(sigma, settings.scatterCluster.minSigma);

                s.variance.at(origIdx) = sigma * sigma;

                std::normal_distribution<double> normDist(0., sigma);
                s.obsPhot.at(origIdx) += normDist(gen);
            }
        }
    }

    // Determine output file
    string outputFile;
    if (!settings.scatterCluster.outputFile.empty())
    {
        outputFile = settings.scatterCluster.outputFile;
    }
    else
    {
        outputFile = settings.global.files.output + ".scatter.out";
    }

    std::ofstream fout(outputFile);

    if (!fout)
    {
        cerr << "***Error: Could not create output file " << outputFile << "***" << endl;
        cerr << "[Exiting...]" << endl;
        exit(-1);
    }

    // Output photometry header
    fout << setprecision(3) << fixed;

    fout << setw(4) << "id" << ' ';

    for (const auto &f : outputFilters)
        fout << setw(6) << f << ' ';

    for (const auto &f : outputFilters)
        fout << setw(6) << ("sig" + f) << ' ';

    fout << setw(7) <<     "mass1" << ' '
         << setw(9) << "massRatio" << ' '
         << setw(5) <<     "stage" << ' '
         << setw(7) <<   "Cmprior" << ' '
         << setw(6) <<    "useDBI" << ' '
         << endl;

    for (size_t i = 0; i < systems.size(); ++i)
    {
        auto &s = systems[i];

        fout << setw(4) << (i + 1) << ' ';

        for (size_t f = 0; f < outputFilters.size(); ++f)
        {
            auto origIdx = filterIndices[f];
            fout << setw(6) << s.obsPhot[origIdx] << ' ';
        }

        for (size_t f = 0; f < outputFilters.size(); ++f)
        {
            auto origIdx = filterIndices[f];
            auto v = s.variance[origIdx];
            fout << setw(6) << (v < 0 ? v : sqrt(v)) << ' ';
        }

        fout << setw(7) << s.primary.mass       << ' '
             << setw(9) << s.getMassRatio()     << ' '
             << setw(5) << s.observedStatus     << ' '
             << setw(7) << s.clustStarPriorDens << ' '
             << setw(6) << 1
             << endl;
    }

    return 0;
}
