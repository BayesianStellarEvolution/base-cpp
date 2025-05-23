#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <vector>

#include "Chain.hpp"

#include "constants.hpp"
#include "densities.hpp"
#include "ifmr.hpp"
#include "IO/SampleWdMass.hpp"
#include "Model.hpp"
#include "Settings.hpp"
#include "WhiteDwarf.hpp"
#include "Utility.hpp"

using std::array;
using std::string;
using std::istringstream;
using std::stringstream;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::unique_ptr;

struct ifmrGridControl
{
    array<double, NPARAMS> priorMean, priorVar;
    double minMag;
    double maxMag;
    int iStart;
    int modelSet;
    array<double, NPARAMS> start; /* starting points for grid evaluations */
    array<double, NPARAMS> end;   /* end points for grid evaluations */
};

/* Used in densities.c. */
array<double, NPARAMS> priorMean, priorVar;

/* TEMPORARY - global variable */
constexpr double const dMass1 = 0.0005;

Settings settings;


/*
 * read control parameters from input stream
 */
static void initIfmrGridControl (Chain *mc, Model &evoModels, struct ifmrGridControl *ctrl, Settings &s)
{
    ctrl->priorMean.at(CARBONICITY) = mc->clust.priorMean.at(CARBONICITY) = s.cluster.priorMeans.carbonicity;
    ctrl->priorVar.at(CARBONICITY) = s.cluster.priorSigma.carbonicity;

    ctrl->priorMean.at(FEH) = s.cluster.priorMeans.Fe_H;
    ctrl->priorVar.at(FEH)  = s.cluster.priorSigma.Fe_H;

    if (ctrl->priorVar.at(FEH) < 0.0)
    {
        ctrl->priorVar.at(FEH) = 0.0;
    }

    ctrl->priorMean.at(MOD) = s.cluster.priorMeans.distMod;
    ctrl->priorVar.at(MOD)  = s.cluster.priorSigma.distMod;
    if (ctrl->priorVar.at(MOD) < 0.0)
    {
        ctrl->priorVar.at(MOD) = 0.0;
    }

    ctrl->priorMean.at(ABS) = s.cluster.priorMeans.Av;
    ctrl->priorVar.at(ABS)  = s.cluster.priorSigma.Av;
    if (ctrl->priorVar.at(ABS) < 0.0)
    {
        ctrl->priorVar.at(ABS) = 0.0;
    }

    ctrl->priorMean.at(AGE) = s.cluster.priorMeans.logAge;
    ctrl->priorVar.at(AGE)  = s.cluster.priorSigma.logAge;

    ctrl->priorMean.at(YYY) = s.cluster.priorMeans.Y;
    ctrl->priorVar.at(YYY)  = s.cluster.priorSigma.Y;
    if (ctrl->priorVar.at(YYY) < 0.0)
    {
        ctrl->priorVar.at(YYY) = 0.0;
    }

    ctrl->priorMean.at(CARBONICITY) = s.cluster.priorMeans.carbonicity;
    ctrl->priorVar.at(CARBONICITY)  = s.cluster.priorSigma.carbonicity;
    if (ctrl->priorVar.at(CARBONICITY) < 0.0)
    {
        ctrl->priorVar.at(CARBONICITY) = 0.0;
    }

    ctrl->priorVar.at(IFMR_INTERCEPT) = 1.0;
    ctrl->priorVar.at(IFMR_SLOPE) = 1.0;

    if (evoModels.IFMR >= 9 && evoModels.IFMR < 12)
        ctrl->priorVar.at(IFMR_QUADCOEF) = 1.0;
    else
        ctrl->priorVar.at(IFMR_QUADCOEF) = 0.0;

    for (auto &var : ctrl->priorVar)
    {
        if (var < 0.0)
        {
            var = 0.0;
        }
        else
        {
            var = var * var;
        }
    }

    // copy values to global variables
    priorVar.at(AGE) = ctrl->priorVar.at(AGE);
    priorVar.at(FEH) = ctrl->priorVar.at(FEH);
    priorVar.at(MOD) = ctrl->priorVar.at(MOD);
    priorVar.at(ABS) = ctrl->priorVar.at(ABS);
    priorVar.at(YYY) = ctrl->priorVar.at(YYY);
    priorVar.at(CARBONICITY) = ctrl->priorVar.at(CARBONICITY);
    priorVar.at(IFMR_INTERCEPT) = ctrl->priorVar.at(IFMR_INTERCEPT);
    priorVar.at(IFMR_SLOPE) = ctrl->priorVar.at(IFMR_SLOPE);
    priorVar.at(IFMR_QUADCOEF) = ctrl->priorVar.at(IFMR_QUADCOEF);

    priorMean.at(FEH) = ctrl->priorMean.at(FEH);
    priorMean.at(MOD) = ctrl->priorMean.at(MOD);
    priorMean.at(ABS) = ctrl->priorMean.at(ABS);
    priorMean.at(YYY) = ctrl->priorMean.at(YYY);
    priorMean.at(CARBONICITY) = ctrl->priorMean.at(CARBONICITY);

    /* prior values for linear IFMR */
    ctrl->priorMean.at(IFMR_SLOPE) = 0.08;
    ctrl->priorMean.at(IFMR_INTERCEPT) = 0.65;
    ctrl->priorMean.at(IFMR_QUADCOEF) = 0.0;
    priorMean.at(IFMR_SLOPE) = ctrl->priorMean.at(IFMR_SLOPE);
    priorMean.at(IFMR_INTERCEPT) = ctrl->priorMean.at(IFMR_INTERCEPT);
    priorMean.at(IFMR_QUADCOEF) = ctrl->priorMean.at(IFMR_QUADCOEF);

    /* open files for reading (data) and writing */

    ctrl->minMag = s.cluster.minMag;
    ctrl->maxMag = s.cluster.maxMag;

    ctrl->iStart = 0;

} // initIfmrGridControl


/*
 * Initialize chain
 */
static void initChain (Chain *mc, const struct ifmrGridControl *ctrl)
{
    int p;

    for (p = 0; p < NPARAMS; p++)
    {
        mc->acceptClust.at(p) = mc->rejectClust.at(p) = 0;
    }

    // If there is no beta in file, initialize everything to prior means
    mc->clust.feh = settings.cluster.starting.Fe_H;
    mc->clust.mod = settings.cluster.starting.Av;
    mc->clust.abs = settings.cluster.starting.Av;
    mc->clust.yyy = settings.cluster.starting.Y;
    mc->clust.age = settings.cluster.starting.logAge;
    mc->clust.carbonicity = settings.cluster.starting.carbonicity;
    mc->clust.ifmrIntercept = ctrl->priorMean.at(IFMR_INTERCEPT);
    mc->clust.ifmrSlope = ctrl->priorMean.at(IFMR_SLOPE);
    mc->clust.ifmrQuadCoef = ctrl->priorMean.at(IFMR_QUADCOEF);

    mc->clust.mean.at(AGE) = ctrl->priorMean.at(AGE);
    mc->clust.mean.at(YYY) = ctrl->priorMean.at(YYY);
    mc->clust.mean.at(MOD) = ctrl->priorMean.at(MOD);
    mc->clust.mean.at(FEH) = ctrl->priorMean.at(FEH);
    mc->clust.mean.at(ABS) = ctrl->priorMean.at(ABS);
    mc->clust.mean.at(CARBONICITY) = ctrl->priorMean.at(CARBONICITY);
    mc->clust.mean.at(IFMR_INTERCEPT) = ctrl->priorMean.at(IFMR_INTERCEPT);
    mc->clust.mean.at(IFMR_SLOPE) = ctrl->priorMean.at(IFMR_SLOPE);
    mc->clust.mean.at(IFMR_QUADCOEF) = ctrl->priorMean.at(IFMR_QUADCOEF);

    for (auto star : mc->stars)
    {
        star.clustStarProposalDens = star.clustStarPriorDens;   // Use prior prob of being clus star

        // find photometry for initial values of currentClust and mc->stars
        if (star.observedStatus == WD)
        {
            star.setMassRatio(0.0);
        }
    }
} // initChain

int main (int argc, char *argv[])
{
    Chain mc;
    struct ifmrGridControl ctrl;

    double fsLike;

    settings.loadSettings (argc, argv);

    if (settings.seed == std::numeric_limits<uint32_t>::max())
    {
        srand(std::time(0));
        settings.seed = rand();

        cout << "Seed: " << settings.seed << endl;
    }

    std::mt19937 gen(settings.seed * uint32_t(2654435761)); // Applies Knuth's multiplicative hash for obfuscation (TAOCP Vol. 3)
    {
        const int warmupIter = 10000;

        cout << "Warming up generator..." << flush;

        for (auto i = 0; i < warmupIter; ++i)
        {
            std::generate_canonical<double, 10>(gen);
        }

        cout << " Done." << endl;
        cout << "Generated " << warmupIter << " values." << endl;
    }

    std::unique_ptr<SampleWdMassBackingStore> sampleWdMassStore;

    sampleWdMassStore.reset(new SampleWdMass_FileBackingStore(settings.files.output));


    Model evoModels = makeModel(settings);

    initIfmrGridControl (&mc, evoModels, &ctrl, settings);

    {
        // Read photometry and calculate fsLike
        vector<double> filterPriorMin;
        vector<double> filterPriorMax;

        vector<string> filterNames;

        std::ifstream rData(settings.files.phot);

        if (!rData)
        {
            cerr << "***Error: Photometry file " << settings.files.phot << " was not found.***" << endl;
            cerr << ".at(Exiting...)" << endl;
            exit (1);
        }

        auto ret = base::utility::readPhotometry (rData, filterPriorMin, filterPriorMax, settings);

        filterNames = ret.first;
        mc.stars = ret.second;

        evoModels.restrictFilters(filterNames, settings.allowInvalidModels);

        if (   settings.cluster.index < 0
               || static_cast<size_t>(settings.cluster.index) > filterNames.size())
        {
            cerr << "*** Error: " << settings.cluster.index << " not a valid magnitude index.  Choose 0 through " << filterNames.size() - 1 << " ***"<< endl;
            cerr << "(Exiting...)" << endl;
            exit (1);
        }

        double logFieldStarLikelihood = 0.0;

        for (size_t filt = 0; filt < filterNames.size(); filt++)
        {
            logFieldStarLikelihood -= log (filterPriorMax.at(filt) - filterPriorMin.at(filt));
        }

        fsLike = mc.stars.size() == 1 ? 0.0 : exp (logFieldStarLikelihood);
    }

    mc.clust.setM_wd_up(settings.whiteDwarf.M_wd_up);

    initChain (&mc, &ctrl);

    auto sampledPars = base::utility::readSampledParams (evoModels, settings);
    cout << "sampledPars.at(0).age    = " << sampledPars.at(0).age << endl;
    cout << "sampledPars.at(last).age = " << sampledPars.back().age << endl;

    /* initialize WD logpost array and WD indices */
    double nWDLogPosts = (int) ceil ((mc.clust.getM_wd_up() - 0.15) / dMass1);

    /********** compile results *********/
    /*** now report sampled masses and parameters ***/

    std::vector<double> us;
    us.resize(sampledPars.size() + 1);

    for (auto &u : us)
    {
        u = std::generate_canonical<double, 53>(gen);
    }

    for (size_t m = 0; m < sampledPars.size(); ++m)
    {
        Cluster internalCluster(mc.clust);
        std::vector<SampleWdMassRecord> records;
        std::vector<double> wdLogPost;

        wdLogPost.resize(nWDLogPosts + 1);

        internalCluster.age = sampledPars.at(m).age;
        internalCluster.feh = sampledPars.at(m).feh;
        internalCluster.mod = sampledPars.at(m).distMod;
        internalCluster.abs = sampledPars.at(m).abs;

        if (evoModels.IFMR >= 4 && evoModels.IFMR < 12)
        {
            internalCluster.ifmrIntercept = sampledPars.at(m).ifmrIntercept;
            internalCluster.ifmrSlope = sampledPars.at(m).ifmrSlope;
        }

        if (evoModels.IFMR >= 9 && evoModels.IFMR < 12)
        {
            internalCluster.ifmrQuadCoef = sampledPars.at(m).ifmrQuadCoef;
        }

        /************ sample WD masses for different parameters ************/
        int im;
        double wdPostSum, maxWDLogPost, mass1;
        double postClusterStar;

        unique_ptr<Isochrone> isochrone(evoModels.mainSequenceEvol->deriveIsochrone(internalCluster.feh, internalCluster.yyy, internalCluster.age));

        auto iteration = sampleWdMassStore->nextIteration();

        for (auto star : mc.stars)
        {
            if (star.observedStatus == WD)
            {
                postClusterStar = 0.0;

                im = 0;

                maxWDLogPost = -HUGE_VAL;

                for (mass1 = 0.15; mass1 < internalCluster.getM_wd_up(); mass1 += dMass1)
                {
                    if (!settings.onlyWDs || mass1 > isochrone->agbTipMass())
                    {
                        /* condition on WD being cluster star */
                        star.primary.mass = mass1;
                        star.setMassRatio(0.0);

                        try
                        {
                            wdLogPost.at(im) = star.logPost (internalCluster, evoModels, *isochrone, settings.modIsParallax);
                            postClusterStar += exp (wdLogPost.at(im));

                            if (wdLogPost.at(im) > maxWDLogPost)
                                maxWDLogPost = wdLogPost.at(im);
                        }
                        catch ( WDBoundsError &e )
                        {
                            cerr << e.what() << endl;

                            wdLogPost.at(im) = -HUGE_VAL;
                        }
                    }
                    else
                    {
                        wdLogPost.at(im) = 0;
                    }

                    im++;
                }

                /* compute the normalizing constant */
                wdPostSum = 0.0;
                im = 0;
                for (mass1 = 0.15; mass1 < internalCluster.getM_wd_up(); mass1 += dMass1)
                {
                    if (!settings.onlyWDs || mass1 > isochrone->agbTipMass())
                    {
                        wdPostSum += exp (wdLogPost.at(im) - maxWDLogPost);
                    }

                    im++;
                }

                /* now sample a particular mass */
                double cumSum = 0.0;
                mass1 = 0.15;
                im = 0;
                while (cumSum < us.at(m) && mass1 < internalCluster.getM_wd_up())
                {
                    if (!settings.onlyWDs || mass1 > isochrone->agbTipMass())
                    {
                        cumSum += exp (wdLogPost.at(im) - maxWDLogPost) / wdPostSum;
                    }

                    mass1 += dMass1;
                    im++;
                }
                mass1 -= dMass1;        /* maybe not necessary */

                postClusterStar *= (internalCluster.getM_wd_up() - 0.15);

                double membership = star.clustStarPriorDens * postClusterStar / (star.clustStarPriorDens * postClusterStar + (1.0 - star.clustStarPriorDens) * fsLike);

                {
                    auto precLogAge = evoModels.mainSequenceEvol->wdPrecLogAge(internalCluster.feh, mass1, internalCluster.yyy);

                    double thisWDMass = intlFinalMassReln (internalCluster, evoModels, mass1);
                    double coolingAge = log10(exp10(internalCluster.age) - exp10(precLogAge));

                    auto teffRadiusPair = evoModels.WDcooling->wdMassToTeffAndRadius (internalCluster.age, internalCluster.carbonicity, precLogAge, thisWDMass);
                    double logTeff = teffRadiusPair.first;

                    auto logg = evoModels.WDAtmosphere->teffToLogg (logTeff, thisWDMass, star.primary.wdType);

                    records.push_back({iteration, star.id,
                                       mass1, membership, precLogAge, coolingAge, logTeff, logg});
                }
            }
        }

        sampleWdMassStore->save(records);
    }

    cout << "Part 2 completed successfully" << endl;

    return 0;
}
