/* Author: Elliot Robinson
 *
 * Outputs CMDs for every line in a .res output file
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <sstream>
#include <vector>

#include "Cluster.hpp"
#include "constants.hpp"
#include "Isochrone.hpp"
#include "Model.hpp"
#include "Settings.hpp"
#include "Star.hpp"
#include "Utility.hpp"

using std::cerr;
using std::endl;
using std::fabs;
using std::stringstream;
using std::ofstream;
using std::string;
using std::unique_ptr;
using std::vector;

const double MIN_MASS1 = 0.15;

struct clustPar
{
    clustPar(double age, double y, double feh, double modulus, double absorption, double carbonicity, double ifmrIntercept, double ifmrSlope, double ifmrQuadCoef, double logPost)
        : age(age), y(y), feh(feh), distMod(modulus), abs(absorption), carbonicity(carbonicity), ifmrIntercept(ifmrIntercept), ifmrSlope(ifmrSlope), ifmrQuadCoef(ifmrQuadCoef), logPost(logPost)
    {}

    double age;
    double y;
    double feh;
    double distMod;
    double abs;
    double carbonicity;
    double ifmrIntercept;
    double ifmrSlope;
    double ifmrQuadCoef;
    double logPost;
};

static vector<clustPar> readSampledParams (Model &evoModels, const GlobalSettings &s)
{
    string line;

    vector<clustPar> sampledPars;

    std::ifstream parsFile;
    parsFile.open(s.files.output + ".res");

    bool hasY, hasCarbonicity;

    getline(parsFile, line); // Parse header

    {
        string sin;
        stringstream in(line);

        in >> sin  // logAge
           >> sin; // Y?

        if (sin == "Y")
        {
            hasY = true;

            in >> sin  // FeH
               >> sin  // Abs
               >> sin  // DistMod
               >> sin; // Carbonicity?

            if (sin == "carbonicity")
                hasCarbonicity = true;
            else
                hasCarbonicity = false;
        }
        else
        {
            hasY = false;

            // This one skips FeH (because it was already read instead of Y)
            in >> sin  // Abs
               >> sin  // DistMod
               >> sin; // Carbonicity?

            if (sin == "carbonicity")
                hasCarbonicity = true;
            else
                hasCarbonicity = false;
        }
    }

    while (getline(parsFile, line))
    {
        stringstream in(line);

        double newAge, newY, newFeh, newMod, newAbs, newCarbonicity, newIInter, newISlope, newIQuad, newLogPost;
        int stage;

        in >> newAge;

        if (hasY)
            in >> newY;
        else
            newY = s.cluster.starting.Y;

        in >> newFeh
           >> newMod
           >> newAbs;

        if (hasCarbonicity)
            in >> newCarbonicity;
        else
            newCarbonicity = s.cluster.starting.carbonicity;

        if (evoModels.IFMR >= 4 && evoModels.IFMR < 12)
        {
            in >> newIInter
               >> newISlope;
        }

        if (evoModels.IFMR >= 9 && evoModels.IFMR < 12)
        {
            in >> newIQuad;
        }

        in >> newLogPost;

        in >> stage;

        if (stage == 3)
            sampledPars.emplace_back(newAge, newY, newFeh, newMod, newAbs, newCarbonicity, newIInter, newISlope, newIQuad, newLogPost);
    }

    parsFile.close();

    return sampledPars;
}


class Application
{
  public:
    Application(Settings s)
        :  evoModels(makeModel(s.global)), settings(s)
    {
        N_WD_MASS1 = GlobalSettings::calcWdMassSteps(s.wdInterpolationPower);
    }

    ~Application()
    {}

    void run();
    Isochrone interpolateIsochrone(const Cluster&, const Isochrone&);
    vector<EvolutionaryPoint> interpolateWDIsochrone(const Cluster&, const Isochrone&, WdAtmosphere);

  private:
    Model evoModels;
    const Settings settings;

    size_t N_WD_MASS1;

    bool hasMSStars = false;
    bool hasWDStars = false;
    bool hasDA = false;
    bool hasDB = false;

    double wdGridMass(int point) const;
};


double Application::wdGridMass(int point) const
{
    double dWdMass1 = (settings.global.whiteDwarf.M_wd_up - MIN_MASS1) / static_cast<double>(N_WD_MASS1);
    return MIN_MASS1 + (dWdMass1 * point);
}


void Application::run()
{
    vector<string> filterNames;

    // for (size_t f = 0; f < filters.size(); ++f)
    // {
    //     try
    //     {
    //         Filters::absCoeffs.at(filters.at(f));
    //     }
    //     catch(std::out_of_range &e)
    //     {
    //         filters.erase(filters.begin() + f);

    //         if (f > 0)
    //             f -= 1;
    //     }
    // }

    {
        // Read photometry and calculate fsLike
        vector<double> filterPriorMin;
        vector<double> filterPriorMax;

        // open files for reading (data) and writing
        // rData implcitly relies on going out of scope to close the photometry file
        // This is awful, but pretty (since this code is, at time of writing, in restricted, anonymous scope
        std::ifstream rData(settings.global.files.phot);

        if (!rData)
        {
            cerr << "***Error: Photometry file " << settings.global.files.phot << " was not found.***" << endl;
            cerr << "[Exiting...]" << endl;
            exit (-1);
        }

        auto ret = base::utility::readPhotometry (rData, filterPriorMin, filterPriorMax, settings.global);
        filterNames = ret.first;

        for (const auto &system : ret.second)
        {
            if (system.observedStatus == MSRG)
            {
                hasMSStars = true;
            }
            else if (system.observedStatus == WD)
            {
                hasWDStars = true;
                if (system.primary.wdType == WdAtmosphere::DA)
                    hasDA = true;
                else if (system.primary.wdType == WdAtmosphere::DB)
                    hasDB = true;
            }
        }

        evoModels.restrictFilters(filterNames, false);
    }

    // evoModels.restrictFilters(filters);
    vector<string> filters = evoModels.mainSequenceEvol->getAvailableFilters();

    Cluster clust;

    auto sampledParams = readSampledParams (evoModels, settings.global);

    double sampledMin = std::numeric_limits<double>::max();
    double sampledMax = std::numeric_limits<double>::lowest();

    for (auto s : sampledParams)
    {
        sampledMin = std::min(sampledMin, exp(s.logPost));
        sampledMax = std::max(sampledMax, exp(s.logPost));
    }

    auto alpha = [=](double x){
        double delta = sampledMax - sampledMin;
        return (exp(x) - sampledMin) / delta;
    };


    string filename = settings.global.files.output + ".sampleCmds";

    ofstream fout(filename);

    if (!fout)
    {
        cerr << "*** Error: Output file " << filename << " could not be opened.***" << endl;
        cerr << "*** [Exiting...]" << endl;
        exit (1);
    }

    // Print file header
    {
        for (auto f : filters)
        {
            fout << base::utility::format << f << ' ';
        }

        fout << base::utility::format << "alpha";

        fout << '\n';
    }

    bool outputMS = !settings.onlyWDs && hasMSStars;
    bool outputWD = !settings.simCluster.noWDs && hasWDStars;

    for (auto s : sampledParams)
    {
        // Standard cluster parameters
        clust.feh = clust.priorMean[FEH] = s.feh;
        clust.mod = clust.priorMean[MOD] = s.distMod;
        clust.abs = clust.priorMean[ABS] = fabs(s.abs);
        clust.age = clust.priorMean[AGE] = s.age;

        // IFMR parameters
        clust.ifmrSlope = clust.priorMean[IFMR_SLOPE] = 0.08;
        clust.ifmrIntercept = clust.priorMean[IFMR_INTERCEPT] = 0.65;

        // FilterSet-conditional parameters
        clust.carbonicity = clust.priorMean[CARBONICITY] = s.carbonicity;
        clust.yyy = clust.priorMean[YYY] = s.y;

        if (evoModels.IFMR <= 10 && evoModels.IFMR < 12)
            clust.ifmrQuadCoef = clust.priorMean[IFMR_QUADCOEF] = 0.0001;
        else
            clust.ifmrQuadCoef = clust.priorMean[IFMR_QUADCOEF] = 0.08;

        clust.setM_wd_up(settings.global.whiteDwarf.M_wd_up);

        // Set Cluster's AGBT_zmass based on the model's agbTipMass
        unique_ptr<Isochrone> isochrone(evoModels.mainSequenceEvol->deriveIsochrone(clust.feh, clust.yyy, clust.age));

        if (outputMS)
        {
            // Interpolate a new isochrone
            auto interpIsochrone = interpolateIsochrone(clust, *isochrone);

            for (auto eep : interpIsochrone.eeps)
            {
                for (auto mag : eep.mags)
                {
                    fout << base::utility::format << mag << ' ';
                }

                fout << base::utility::format << alpha(s.logPost);

                fout << '\n';
            }
        }

        if (outputWD)
        {
            vector<EvolutionaryPoint> wdEeps;

            if (hasDA)
            {
                auto daEeps = interpolateWDIsochrone(clust, *isochrone, WdAtmosphere::DA);
                wdEeps.insert(wdEeps.end(), daEeps.begin(), daEeps.end());
            }

            if (hasDB)
            {
                auto dbEeps = interpolateWDIsochrone(clust, *isochrone, WdAtmosphere::DB);
                wdEeps.insert(wdEeps.end(), dbEeps.begin(), dbEeps.end());
            }

            for (auto eep : wdEeps)
            {
                for (auto mag : eep.mags)
                {
                    fout << base::utility::format << mag << ' ';
                }

                fout << base::utility::format << alpha(s.logPost);

                fout << '\n';
            }
        }
    }

    fout.close();
}


Isochrone Application::interpolateIsochrone(const Cluster &clust, const Isochrone &isochrone)
{
    vector<EvolutionaryPoint> eeps;

    if (settings.interpolateCMD)
    {
        const int isoIncrem =
            settings.eepInterpolationPower == 0 ? 80
                                                : settings.eepInterpolationPower > 0 ? 32 << settings.eepInterpolationPower
                                                                                     : 32 >> -settings.eepInterpolationPower;

        auto isoSize = isochrone.eeps.size();

        for (size_t m = 0; m < isoSize - 1; ++m)
        {
            double dIsoMass = isochrone.eeps[m + 1].mass - isochrone.eeps[m].mass;

            if (dIsoMass < 0.0) continue;

            if (isochrone.eeps[m].mass > isochrone.agbTipMass())
                continue;

            double dMass = dIsoMass / isoIncrem;

            for (auto k = 0; k < isoIncrem; ++k)
            {
                double mass = isochrone.eeps[m].mass + (k * dMass);

                if (mass > isochrone.agbTipMass())
                    break;

                StellarSystem star;
                star.primary.mass   = mass;
                star.secondary.mass = 0.0;

                eeps.emplace_back(m * isoIncrem + k, star.primary.mass, star.deriveCombinedMags(clust, evoModels, isochrone, settings.global.modIsParallax));
            }
        }
    }
    else
    {
        for (size_t e = 0; e < isochrone.eeps.size(); ++e)
        {
            double mass = isochrone.eeps.at(e).mass;

            if (mass > isochrone.agbTipMass())
                continue;

            StellarSystem star;
            star.primary.mass   = mass;
            star.secondary.mass = 0.0;

            eeps.emplace_back(e, star.primary.mass, star.deriveCombinedMags(clust, evoModels, isochrone, settings.global.modIsParallax));
        }
    }

    return {clust.age, eeps};
}


vector<EvolutionaryPoint> Application::interpolateWDIsochrone(const Cluster &clust, const Isochrone &isochrone, WdAtmosphere wdType)
{
    vector<EvolutionaryPoint> eeps;

    if (settings.interpolateCMD)
    {
        for (size_t j = 0; j < N_WD_MASS1; ++j)
        {
            double mass = wdGridMass(j);

            if (mass <= isochrone.agbTipMass() || mass > settings.global.whiteDwarf.M_wd_up)
                continue;

            {
                StellarSystem star;
                star.primary.mass   = mass;
                star.primary.wdType = wdType;
                star.secondary.mass = 0.0;

                eeps.emplace_back(j, star.primary.mass, star.deriveCombinedMags(clust, evoModels, isochrone, settings.global.modIsParallax));
            }
        }
    }
    else
    {
        int e = 0;
        double mass = isochrone.agbTipMass();

        while (mass < settings.global.whiteDwarf.M_wd_up)
        {
            mass = isochrone.agbTipMass() + ++e * 0.05;

            if (mass > settings.global.whiteDwarf.M_wd_up)
                break;

            StellarSystem star;
            star.primary.mass   = mass;
            star.primary.wdType = wdType;
            star.secondary.mass = 0.0;

            eeps.emplace_back(e, star.primary.mass, star.deriveCombinedMags(clust, evoModels, isochrone, settings.global.modIsParallax));
        }
    }

    return eeps;
}


int main (int argc, char *argv[])
{
    Settings settings;

    settings.loadSettings (argc, argv);

    if (settings.global.seed == std::numeric_limits<uint32_t>::max())
    {
        srand(std::time(0));
        settings.global.seed = rand();

        cout << "Seed: " << settings.global.seed << endl;
    }

    Application(settings).run();

    return 0;
}
