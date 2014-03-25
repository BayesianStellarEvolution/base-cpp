#include <array>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#include <boost/format.hpp>

#include "Chain.hpp"

#include "constants.hpp"
#include "densities.hpp"
#include "Model.hpp"
#include "Settings.hpp"
#include "WhiteDwarf.hpp"
#include "Utility.hpp"

using std::array;
using std::string;
using std::istringstream;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::flush;

const int N_AGE = 30;
const int N_FEH = 1;
const int N_MOD = 1;
const int N_ABS = 1;
const int N_Y = 1;
const int N_IFMR_INT = 10;
const int N_IFMR_SLOPE = 10;
const int N_GRID = (N_AGE * N_FEH * N_MOD * N_ABS * N_Y * N_IFMR_INT * N_IFMR_SLOPE);
const int MASTER = 0;       /* taskid of first process */

const int ALLOC_CHUNK = 1;

struct ifmrGridControl
{
    double initialAge;
    array<double, NPARAMS> priorMean, priorVar;
    double minMag;
    double maxMag;
    int iMag;
    int iStart;
    int modelSet;
    array<double, FILTS> filterPriorMin, filterPriorMax;
    int numFilts;
    int nSamples;
    array<double, NPARAMS> start; /* starting points for grid evaluations */
    array<double, NPARAMS> end;   /* end points for grid evaluations */
};

/* For posterior evaluation on a grid */
struct clustPar
{
    clustPar(double age, double feh, double modulus, double absorption, double ifmrIntercept, double ifmrSlope, double ifmrQuadCoef)
        : age(age), FeH(feh), modulus(modulus), absorption(absorption), ifmrIntercept(ifmrIntercept), ifmrSlope(ifmrSlope), ifmrQuadCoef(ifmrQuadCoef)
    {}

    double age;
    double FeH;
    double modulus;
    double absorption;
    double ifmrIntercept;
    double ifmrSlope;
    double ifmrQuadCoef;
};

typedef struct
{
    array<double, FILTS> obsPhot, variance;
    double clustStarPriorDens;  /* cluster membership prior probability */
} obsStar;

/* declare global variables */
array<double, FILTS> filterPriorMin;
array<double, FILTS> filterPriorMax;

/* Used in densities.c. */
array<double, NPARAMS> priorMean, priorVar;

/* Used by a bunch of different functions. */
vector<int> filters;

/* TEMPORARY - global variable */
constexpr double const dMass1 = 0.0005;

Settings settings;


/*
 * read control parameters from input stream
 */
static void initIfmrGridControl (Chain *mc, Model &evoModels, struct ifmrGridControl *ctrl, Settings &s)
{
    ctrl->numFilts = 0;

    if (s.whiteDwarf.wdModel == WdModel::MONTGOMERY)
    {
        ctrl->priorMean.at(CARBONICITY) = mc->clust.carbonicity = mc->clust.priorMean.at(CARBONICITY) = settings.cluster.carbonicity;
        ctrl->priorVar.at(CARBONICITY) = settings.cluster.sigma.carbonicity;
    }
    else
    {
        ctrl->priorMean.at(CARBONICITY) = mc->clust.carbonicity = mc->clust.priorMean.at(CARBONICITY) = 0.0;
        ctrl->priorVar.at(CARBONICITY) = 0.0;
    }


    ctrl->priorMean.at(FEH) = settings.cluster.Fe_H;
    ctrl->priorVar.at(FEH) = settings.cluster.sigma.Fe_H;
    if (ctrl->priorVar.at(FEH) < 0.0)
    {
        ctrl->priorVar.at(FEH) = 0.0;
    }

    ctrl->priorMean.at(MOD) = settings.cluster.distMod;
    ctrl->priorVar.at(MOD) = settings.cluster.sigma.distMod;
    if (ctrl->priorVar.at(MOD) < 0.0)
    {
        ctrl->priorVar.at(MOD) = 0.0;
    }

    ctrl->priorMean.at(ABS) = settings.cluster.Av;
    ctrl->priorVar.at(ABS) = settings.cluster.sigma.Av;
    if (ctrl->priorVar.at(ABS) < 0.0)
    {
        ctrl->priorVar.at(ABS) = 0.0;
    }

    ctrl->initialAge = settings.cluster.logClusAge;
    ctrl->priorVar.at(AGE) = 1.0;

    ctrl->priorVar.at(IFMR_INTERCEPT) = 1.0;
    ctrl->priorVar.at(IFMR_SLOPE) = 1.0;
    if (evoModels.IFMR >= 9)
        ctrl->priorVar.at(IFMR_QUADCOEF) = 1.0;
    else
        ctrl->priorVar.at(IFMR_QUADCOEF) = 0.0;

    // copy values to global variables
    priorVar.at(AGE) = ctrl->priorVar.at(AGE);
    priorVar.at(FEH) = ctrl->priorVar.at(FEH);
    priorVar.at(MOD) = ctrl->priorVar.at(MOD);
    priorVar.at(ABS) = ctrl->priorVar.at(ABS);
    priorVar.at(IFMR_INTERCEPT) = ctrl->priorVar.at(IFMR_INTERCEPT);
    priorVar.at(IFMR_SLOPE) = ctrl->priorVar.at(IFMR_SLOPE);
    priorVar.at(IFMR_QUADCOEF) = ctrl->priorVar.at(IFMR_QUADCOEF);

    priorMean.at(FEH) = ctrl->priorMean.at(FEH);
    priorMean.at(MOD) = ctrl->priorMean.at(MOD);
    priorMean.at(ABS) = ctrl->priorMean.at(ABS);

    /* prior values for linear IFMR */
    ctrl->priorMean.at(IFMR_SLOPE) = 0.08;
    ctrl->priorMean.at(IFMR_INTERCEPT) = 0.65;
    ctrl->priorMean.at(IFMR_QUADCOEF) = 0.0;
    priorMean.at(IFMR_SLOPE) = ctrl->priorMean.at(IFMR_SLOPE);
    priorMean.at(IFMR_INTERCEPT) = ctrl->priorMean.at(IFMR_INTERCEPT);
    priorMean.at(IFMR_QUADCOEF) = ctrl->priorMean.at(IFMR_QUADCOEF);

    /* open model file, choose model set, and load models */

    if (s.mainSequence.msRgbModel == MsModel::CHABHELIUM)
    {
        scanf ("%lf %lf", &ctrl->priorMean.at(YYY), &ctrl->priorVar.at(YYY));

        if (ctrl->priorVar.at(YYY) < 0.0)
        {
            ctrl->priorVar.at(YYY) = 0.0;
        }
    }
    else
    {
        ctrl->priorMean.at(YYY) = 0.0;
        ctrl->priorVar.at(YYY) = 0.0;
    }
    priorVar.at(YYY) = ctrl->priorVar.at(YYY);
    priorMean.at(YYY) = ctrl->priorMean.at(YYY);

    /* open files for reading (data) and writing */

    ctrl->minMag = settings.cluster.minMag;
    ctrl->maxMag = settings.cluster.maxMag;
    ctrl->iMag = settings.cluster.index;
    if (ctrl->iMag < 0 || ctrl->iMag > FILTS)
    {
        cerr << "***Error: " << ctrl->iMag << " not a valid magnitude index.  Choose 0, 1,or 2.***" << endl;
        cerr << ".at(Exiting...)" << endl;
        exit (1);
    }

    ctrl->iStart = 0;

    /* Initialize filter prior mins and maxes */
    int j;

    for (j = 0; j < FILTS; j++)
    {
        ctrl->filterPriorMin.at(j) = 1000;
        ctrl->filterPriorMax.at(j) = -1000;
    }
} // initIfmrGridControl

void readCmdData (Chain &mc, struct ifmrGridControl &ctrl, const Model &evoModels)
{
    string line, pch;

    std::ifstream rData;
    rData.open(settings.files.phot);
    if (!rData)
    {
        cerr << "***Error: Photometry file " << settings.files.phot << " was not found.***" << endl;
        cerr << ".at(Exiting...)" << endl;
        exit (1);
    }

    //Parse the header of the file to determine which filters are being used
    getline(rData, line);  // Read in the header line

    istringstream header(line); // Ignore the first token (which is "id")

    header >> pch;

    while (!header.eof())
    {
        header >> pch;

        if (pch == "sig")
            break;                      // and check to see if they are 'sig'.  If they are, there are no more filters

        for (int filt = 0; filt < FILTS; filt++)
        {                               // Otherwise check to see what this filter's name is
            if (pch == evoModels.filterSet->getFilterName(filt))
            {
                ctrl.numFilts++;
                filters.push_back(filt);
                const_cast<Model&>(evoModels).numFilts++;
                break;
            }
        }
    }

    // This loop reads in photometry data
    // It also reads a best guess for the mass
    mc.stars.clear();

    while (!rData.eof())
    {
        getline(rData, line);

        if (rData.eof())
            break;

        mc.stars.emplace_back(line, ctrl.numFilts);

        for (int i = 0; i < ctrl.numFilts; i++)
        {
            if (mc.stars.back().obsPhot.at(i) < ctrl.filterPriorMin.at(i))
            {
                ctrl.filterPriorMin.at(i) = mc.stars.back().obsPhot.at(i);
            }

            if (mc.stars.back().obsPhot.at(i) > ctrl.filterPriorMax.at(i))
            {
                ctrl.filterPriorMax.at(i) = mc.stars.back().obsPhot.at(i);
            }

            filterPriorMin.at(i) = ctrl.filterPriorMin.at(i);
            filterPriorMax.at(i) = ctrl.filterPriorMax.at(i);
        }

        if (!(mc.stars.back().observedStatus == WD || (mc.stars.back().obsPhot.at(ctrl.iMag) >= ctrl.minMag && mc.stars.back().obsPhot.at(ctrl.iMag) <= ctrl.maxMag)))
        {
            mc.stars.pop_back();
        }
    }

    rData.close();
} /* readCmdData */

/*
 * Read sampled params
 */
static void readSampledParams (struct ifmrGridControl *ctrl, vector<clustPar> &sampledPars, Model &evoModels)
{
    string line;
    std::ifstream parsFile;
    parsFile.open(settings.files.output + ".res");

    getline(parsFile, line); // Eat the header

    while (!parsFile.eof())
    {
        double newAge, newFeh, newMod, newAbs, newIInter, newISlope, newIQuad, ignore;
        newAge = newFeh = newMod = newAbs = newIInter = newISlope = newIQuad = 0.0;

        parsFile >> newAge
                 >> newFeh
                 >> newMod
                 >> newAbs;

        if (evoModels.IFMR >= 4)
        {
            parsFile >> newIInter
                     >> newISlope;
        }

        if (evoModels.IFMR >= 9)
        {
            parsFile >> newIQuad;
        }

        parsFile >> ignore; // logPost

        if (!parsFile.eof())
        {
            sampledPars.emplace_back(newAge, newFeh, newMod, newAbs, newIInter, newISlope, newIQuad);
        }
    }

    parsFile.close();

    ctrl->nSamples = sampledPars.size();
}


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
    mc->clust.feh = ctrl->priorMean.at(FEH);
    mc->clust.mod = ctrl->priorMean.at(MOD);
    mc->clust.abs = ctrl->priorMean.at(ABS);
    mc->clust.yyy = ctrl->priorMean.at(YYY);
    mc->clust.age = ctrl->initialAge;
    mc->clust.ifmrIntercept = ctrl->priorMean.at(IFMR_INTERCEPT);
    mc->clust.ifmrSlope = ctrl->priorMean.at(IFMR_SLOPE);
    mc->clust.ifmrQuadCoef = ctrl->priorMean.at(IFMR_QUADCOEF);
    mc->clust.mean.at(AGE) = ctrl->initialAge;
    mc->clust.mean.at(YYY) = ctrl->priorMean.at(YYY);
    mc->clust.mean.at(MOD) = ctrl->priorMean.at(MOD);
    mc->clust.mean.at(FEH) = ctrl->priorMean.at(FEH);
    mc->clust.mean.at(ABS) = ctrl->priorMean.at(ABS);
    mc->clust.mean.at(IFMR_INTERCEPT) = ctrl->priorMean.at(IFMR_INTERCEPT);
    mc->clust.mean.at(IFMR_SLOPE) = ctrl->priorMean.at(IFMR_SLOPE);
    mc->clust.mean.at(IFMR_QUADCOEF) = ctrl->priorMean.at(IFMR_QUADCOEF);

    for (auto star : mc->stars)
    {
        star.clustStarProposalDens = star.clustStarPriorDens;   // Use prior prob of being clus star

        star.primary.wdType = WdAtmosphere::DA;
        star.secondary.wdType = WdAtmosphere::DA;

        // find photometry for initial values of currentClust and mc->stars
        if (star.observedStatus == WD)
        {
            star.setMassRatio(0.0);
        }
    }
} // initChain

int main (int argc, char *argv[])
{
    int filt, nWDs = 0;

    Chain mc;
    struct ifmrGridControl ctrl;

    double fsLike;

    vector<clustPar> sampledPars;

    settings.fromCLI (argc, argv);
    if (!settings.files.config.empty())
    {
        settings.fromYaml (settings.files.config);
    }
    else
    {
        settings.fromYaml ("base9.yaml");
    }

    settings.fromCLI (argc, argv);

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

    Model evoModels = makeModel(settings);

    initIfmrGridControl (&mc, evoModels, &ctrl, settings);

    readCmdData (mc, ctrl, evoModels);

    mc.clust.setM_wd_up(settings.whiteDwarf.M_wd_up);

    vector<obsStar> obs;
    obs.resize(mc.stars.size());

    vector<int> starStatus;
    starStatus.reserve(mc.stars.size());

    for (decltype(mc.stars.size()) i = 0; i < mc.stars.size(); i++)
    {
        for (filt = 0; filt < ctrl.numFilts; filt++)
        {
            obs.at(i).obsPhot.at(filt) = mc.stars.at(i).obsPhot.at(filt);
            obs.at(i).variance.at(filt) = mc.stars.at(i).variance.at(filt);
        }
        obs.at(i).clustStarPriorDens = mc.stars.at(i).clustStarPriorDens;
        starStatus.at(i) = mc.stars.at(i).observedStatus;

        if (starStatus.at(i) == WD)
        {
            nWDs++;
        }
    }

    evoModels.numFilts = ctrl.numFilts;

    initChain (&mc, &ctrl);

    double logFieldStarLikelihood = 0.0;

    for (filt = 0; filt < ctrl.numFilts; filt++)
    {
        logFieldStarLikelihood -= log (ctrl.filterPriorMax.at(filt) - ctrl.filterPriorMin.at(filt));
    }
    fsLike = exp (logFieldStarLikelihood);

    readSampledParams (&ctrl, sampledPars, evoModels);
    cout << "sampledPars.at(0).age    = " << sampledPars.at(0).age << endl;
    cout << "sampledPars.at(last).age = " << sampledPars.back().age << endl;

    /* initialize WD logpost array and WD indices */
    double nWDLogPosts = (int) ceil ((mc.clust.getM_wd_up() - 0.15) / dMass1);

    /********** compile results *********/
    /*** now report sampled masses and parameters ***/

    // Open the file
    string filename = settings.files.output + ".massSamples";

    std::ofstream massSampleFile;
    massSampleFile.open(filename);
    if (!massSampleFile)
    {
        cerr << "***Error: File " << filename << " was not available for writing.***" << endl;
        cerr << ".at(Exiting...)" << endl;
        exit (1);
    }

    filename += ".membership";

    std::ofstream membershipFile;
    membershipFile.open(filename);
    if (!membershipFile)
    {
        cerr << "***Error: File " << filename << " was not available for writing.***" << endl;
        cerr << ".at(Exiting...)" << endl;
        exit (1);
    }
    
    std::vector<double> us;
    us.resize(ctrl.nSamples + 1);

    for (auto &u : us)
    {
        u = std::generate_canonical<double, 53>(gen);
    }

    for (int m = 0; m < ctrl.nSamples; ++m)
    {
        Cluster internalCluster(mc.clust);
        std::vector<double> wdMass, clusMemPost;
        std::vector<double> wdLogPost;

        wdLogPost.resize(nWDLogPosts + 1);

        internalCluster.age = sampledPars.at(m).age;
        internalCluster.feh = sampledPars.at(m).FeH;
        internalCluster.mod = sampledPars.at(m).modulus;
        internalCluster.abs = sampledPars.at(m).absorption;

        if (evoModels.IFMR >= 4)
        {
            internalCluster.ifmrIntercept = sampledPars.at(m).ifmrIntercept;
            internalCluster.ifmrSlope = sampledPars.at(m).ifmrSlope;
        }

        if (evoModels.IFMR >= 9)
        {
            internalCluster.ifmrQuadCoef = sampledPars.at(m).ifmrQuadCoef;
        }

        /************ sample WD masses for different parameters ************/
        int iWD = 0;
        int im;
        double wdPostSum, maxWDLogPost, mass1;
        double postClusterStar;

        internalCluster.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(filters, internalCluster.feh, internalCluster.yyy, internalCluster.age);

        for (auto star : mc.stars)
        {
            if (star.observedStatus == WD)
            {
                postClusterStar = 0.0;

                im = 0;

                for (mass1 = 0.15; mass1 < internalCluster.getM_wd_up(); mass1 += dMass1)
                {
                    /* condition on WD being cluster star */
                    star.primary.mass = mass1;
                    star.setMassRatio(0.0);

                    try
                    {
                        wdLogPost.at(im) = star.logPost (internalCluster, evoModels, filters);
                        postClusterStar += exp (wdLogPost.at(im));
                    }
                    catch ( WDBoundsError &e )
                    {
                        cerr << e.what() << endl;

                        wdLogPost.at(im) = -HUGE_VAL;
                    }

                    im++;
                }
                im = 0;

                /* compute the maximum value */
                maxWDLogPost = wdLogPost.at(0);
                for (mass1 = 0.15; mass1 < internalCluster.getM_wd_up(); mass1 += dMass1)
                {
                    if (wdLogPost.at(im) > maxWDLogPost)
                        maxWDLogPost = wdLogPost.at(im);
                    im++;
                }

                /* compute the normalizing constant */
                wdPostSum = 0.0;
                im = 0;
                for (mass1 = 0.15; mass1 < internalCluster.getM_wd_up(); mass1 += dMass1)
                {
                    wdPostSum += exp (wdLogPost.at(im) - maxWDLogPost);
                    im++;
                }

                /* now sample a particular mass */
                double cumSum = 0.0;
                mass1 = 0.15;
                im = 0;
                while (cumSum < us.at(m) && mass1 < internalCluster.getM_wd_up())
                {
                    cumSum += exp (wdLogPost.at(im) - maxWDLogPost) / wdPostSum;
                    mass1 += dMass1;
                    im++;
                }
                mass1 -= dMass1;        /* maybe not necessary */

                wdMass.push_back(mass1);

                postClusterStar *= (internalCluster.getM_wd_up() - 0.15);

                clusMemPost.push_back(star.clustStarPriorDens * postClusterStar / (star.clustStarPriorDens * postClusterStar + (1.0 - star.clustStarPriorDens) * fsLike));
                iWD++;
            }
        }

        massSampleFile << boost::format("%10.6f") % sampledPars.at(m).age
                       << boost::format("%10.6f") % sampledPars.at(m).FeH
                       << boost::format("%10.6f") % sampledPars.at(m).modulus
                       << boost::format("%10.6f") % sampledPars.at(m).absorption;

        if (evoModels.IFMR >= 4)
        {
            massSampleFile << boost::format("%10.6f") % sampledPars.at(m).ifmrIntercept
                           << boost::format("%10.6f") % sampledPars.at(m).ifmrSlope;
        }

        if (evoModels.IFMR >= 9)
        {
            massSampleFile << boost::format("%10.6f") % sampledPars.at(m).ifmrQuadCoef;
        }

        for (int j = 0; j < nWDs; j++)
        {
            massSampleFile << boost::format("%10.6f") % wdMass.at(j);
            membershipFile << boost::format("%10.6f") % clusMemPost.at(j);
        }

        massSampleFile << endl;
        membershipFile << endl;
    }

    massSampleFile.close();
    membershipFile.close();

    cout << "Part 2 completed successfully" << endl;

    return 0;
}
