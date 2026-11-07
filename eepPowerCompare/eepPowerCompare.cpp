// Copyright (c) Embry-Riddle Aeronautical University. All rights reserved.
//
// Compare log-likelihood values across different EEP interpolation powers
// to assess whether lower powers preserve the likelihood surface shape.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "Cluster.hpp"
#include "marg.hpp"
#include "Model.hpp"
#include "Settings.hpp"
#include "Star.hpp"
#include "Utility.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::vector;

const double TWO_M_PI = 2 * M_PI;

double computeLogLikelihood(
    Cluster& clust,
    Model& evoModels,
    vector<StellarSystem>& systems,
    const double* sysVars,
    const double* sysVar2,
    const double* sysObs,
    const double* sysFluxLo,
    const double* sysFluxHi,
    size_t howManyFiltsAligned,
    size_t howManyFilts,
    base::utility::ThreadPool& pool,
    bool noBinaries,
    int eepPower)
{
    Isochrone* isochrone = evoModels.mainSequenceEvol->deriveIsochrone(clust.feh, clust.yyy, clust.age);

    vector<double> post;
    if (noBinaries) {
        post = margEvolveNoBinaries(clust, evoModels, *isochrone, pool,
                                    sysVars, sysVar2, sysObs,
                                    systems.size(), howManyFiltsAligned, howManyFilts,
                                    false, eepPower);
    } else {
        post = margEvolveWithBinary(clust, evoModels, *isochrone, pool,
                                    systems,
                                    sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi,
                                    howManyFiltsAligned, howManyFilts,
                                    false, eepPower);
    }

    delete isochrone;

    double logLike = 0.0;
    for (size_t i = 0; i < systems.size(); ++i) {
        if (post[i] > 0.0) {
            logLike += std::log(post[i]);
        } else {
            logLike += -1e10;
        }
    }

    return logLike;
}

int main(int argc, char** argv) {
    Settings settings;

    try {
        settings.loadSettings(argc, argv);
    } catch (const std::exception& e) {
        cerr << "Error loading settings: " << e.what() << endl;
        return 1;
    }

    // Load models
    cerr << "Loading models..." << std::flush;
    Model evoModels = makeModel(settings.global);
    cerr << " Done." << endl;

    // Load photometry
    vector<double> filterPriorMin, filterPriorMax;
    std::ifstream photFile(settings.global.files.phot);
    if (!photFile) {
        cerr << "Error: Cannot open photometry file: " << settings.global.files.phot << endl;
        return 1;
    }

    auto photResult = base::utility::readPhotometry(photFile, filterPriorMin, filterPriorMax, settings.global);
    auto filterNames = photResult.first;
    auto systems = photResult.second;
    cerr << "Loaded " << systems.size() << " stars with " << filterNames.size() << " filters" << endl;

    if (systems.empty()) {
        cerr << "Error: No systems loaded" << endl;
        return 1;
    }

    // Restrict filters to those in photometry
    evoModels.restrictFilters(filterNames, false);

    // Set up cluster from starting values
    Cluster baseClust;
    baseClust.feh = settings.global.cluster.starting.Fe_H;
    baseClust.mod = settings.global.cluster.starting.distMod;
    baseClust.abs = settings.global.cluster.starting.Av;
    baseClust.age = settings.global.cluster.starting.logAge;
    baseClust.carbonicity = settings.global.cluster.starting.carbonicity;
    baseClust.yyy = settings.global.cluster.starting.Y;
    baseClust.ifmrSlope = 0.08;
    baseClust.ifmrIntercept = 0.65;
    baseClust.ifmrQuadCoef = 0.0001;
    baseClust.setM_wd_up(settings.global.whiteDwarf.M_wd_up);
    baseClust.varScale = 1.0;

    // Prepare aligned arrays
    size_t howManyFilts = filterNames.size();
    size_t howManyFiltsAligned = ((howManyFilts + 1) & ~0x1);
    size_t totalSize = systems.size() * howManyFiltsAligned;

    vector<double> sysVars(totalSize), sysVar2(totalSize), sysObs(totalSize);
    vector<double> sysFluxLo(totalSize), sysFluxHi(totalSize);

    size_t idx = 0;
    for (const auto& sys : systems) {
        for (size_t k = 0; k < howManyFiltsAligned; ++k, ++idx) {
            if (k < howManyFilts && sys.variance[k] > 1e-10) {
                double var = sys.variance[k] * baseClust.varScale;
                sysVars[idx] = 1.0 / var;
                sysVar2[idx] = std::log(TWO_M_PI * var);
                sysObs[idx] = sys.obsPhot[k];
                double nSDTerm = 3.0 * std::sqrt(var);
                sysFluxLo[idx] = std::pow(10.0, (sys.obsPhot[k] - nSDTerm) / -2.5);
                sysFluxHi[idx] = std::pow(10.0, (sys.obsPhot[k] + nSDTerm) / -2.5);
            } else {
                sysVars[idx] = 0.0;
                sysVar2[idx] = 0.0;
                sysObs[idx] = 0.0;
                sysFluxLo[idx] = 0.0;
                sysFluxHi[idx] = 0.0;
            }
        }
    }

    base::utility::ThreadPool pool(settings.threads);

    // Define parameter ranges (centered on starting values)
    double ageCenter = baseClust.age;
    double ageSpan = 0.2;
    double fehCenter = baseClust.feh;
    double fehSpan = 0.2;
    double modCenter = baseClust.mod;
    double modSpan = 0.4;

    int nPowers = settings.maxPower - settings.minPower + 1;
    int nGridPoints = settings.gridPoints * settings.gridPoints * settings.gridPoints;

    vector<vector<double>> logLikes(nPowers, vector<double>(nGridPoints));

    cerr << "Computing likelihoods on " << nGridPoints << " grid points for powers "
         << settings.minPower << " to " << settings.maxPower << "..." << endl;

    int gridIdx = 0;
    for (int iAge = 0; iAge < settings.gridPoints; ++iAge) {
        double age = ageCenter - ageSpan/2 + ageSpan * iAge / (settings.gridPoints - 1);

        for (int iFeh = 0; iFeh < settings.gridPoints; ++iFeh) {
            double feh = fehCenter - fehSpan/2 + fehSpan * iFeh / (settings.gridPoints - 1);

            for (int iMod = 0; iMod < settings.gridPoints; ++iMod) {
                double mod = modCenter - modSpan/2 + modSpan * iMod / (settings.gridPoints - 1);

                Cluster clust = baseClust;
                clust.age = age;
                clust.feh = feh;
                clust.mod = mod;

                for (int power = settings.minPower; power <= settings.maxPower; ++power) {
                    int powerIdx = power - settings.minPower;
                    logLikes[powerIdx][gridIdx] = computeLogLikelihood(
                        clust, evoModels, systems,
                        sysVars.data(), sysVar2.data(), sysObs.data(),
                        sysFluxLo.data(), sysFluxHi.data(),
                        howManyFiltsAligned, howManyFilts,
                        pool, settings.noBinaries, power);
                }
                ++gridIdx;
            }
        }
        cerr << "  Age slice " << (iAge + 1) << "/" << settings.gridPoints << " complete" << endl;
    }

    // Compute statistics
    int refIdx = settings.referencePower - settings.minPower;

    cout << "\n=== EEP Power Comparison Results ===" << endl;
    cout << "Reference power: " << settings.referencePower << endl;
    cout << "Grid points: " << nGridPoints << " (" << settings.gridPoints << "^3)" << endl;
    cout << "Stars: " << systems.size() << ", Filters: " << howManyFilts << endl;
    cout << "\n";

    cout << std::setw(8) << "Power"
         << std::setw(15) << "Mean Diff"
         << std::setw(15) << "Std Dev"
         << std::setw(15) << "Min Diff"
         << std::setw(15) << "Max Diff"
         << std::setw(12) << "Shape OK?"
         << endl;
    cout << std::string(80, '-') << endl;

    for (int power = settings.minPower; power <= settings.maxPower; ++power) {
        int powerIdx = power - settings.minPower;

        vector<double> diffs(nGridPoints);
        for (int i = 0; i < nGridPoints; ++i) {
            diffs[i] = logLikes[powerIdx][i] - logLikes[refIdx][i];
        }

        double sum = 0.0;
        for (double d : diffs) sum += d;
        double mean = sum / nGridPoints;

        double sumSq = 0.0;
        for (double d : diffs) sumSq += (d - mean) * (d - mean);
        double stdDev = std::sqrt(sumSq / nGridPoints);

        double minDiff = *std::min_element(diffs.begin(), diffs.end());
        double maxDiff = *std::max_element(diffs.begin(), diffs.end());

        bool shapeOK = (stdDev < 0.1);

        cout << std::setw(8) << power
             << std::setw(15) << std::fixed << std::setprecision(4) << mean
             << std::setw(15) << stdDev
             << std::setw(15) << minDiff
             << std::setw(15) << maxDiff
             << std::setw(12) << (shapeOK ? "YES" : "NO")
             << endl;
    }

    cout << "\nInterpretation:" << endl;
    cout << "  Mean Diff: constant offset (OK - cancels in MCMC ratios)" << endl;
    cout << "  Std Dev: shape distortion (BAD - changes posteriors)" << endl;
    cout << "  Shape OK: Std Dev < 0.1" << endl;

    return 0;
}
