#include <functional>
#include <iostream>

#include "Chain.hpp"
#include "Cluster.hpp"
#include "marg.hpp"
#include "mpiMcmc.hpp"
#include "MpiMcmcApplication.hpp"
#include "Settings.hpp"
#include "Star.hpp"
#include "samplers.cpp"
#include "WhiteDwarf.hpp"

using std::array;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::function;
using std::mutex;
using std::tuple;
using std::unique_ptr;
using std::vector;

using namespace std::placeholders;

const double TWO_M_PI = 2 * M_PI;

void ensurePriors(const Settings &s, const DualPopCluster &clust)
{
    // Clust A carbonicity
    if (clust.clust[0].carbonicity < 0.0)
        throw InvalidCluster("Low carbonicity, Clust A");
    else if (clust.clust[0].carbonicity > 1.0)
        throw InvalidCluster("High carbonicity, Clust A");

    // Clust B carbonicity
    else if (clust.clust[1].carbonicity < 0.0)
        throw InvalidCluster("Low carbonicity, Clust B");
    else if (clust.clust[1].carbonicity > 1.0)
        throw InvalidCluster("High carbonicity, Clust B");

    // Clust A Y
    else if (clust.clust[0].yyy < s.multiPopMcmc.YA_lo)
        throw InvalidCluster("Low Y, Clust A");
    else if (clust.clust[0].yyy > s.multiPopMcmc.YA_hi)
        throw InvalidCluster("High Y, Clust A");

    // Clust B Y
    else if (clust.clust[1].yyy < clust.clust[0].yyy)
        throw InvalidCluster("Low Y, Clust B");
    else if (clust.clust[1].yyy > s.multiPopMcmc.YB_hi)
        throw InvalidCluster("High Y, Clust B");

    // Lambda
    else if (clust.lambda < 0.0)
        throw InvalidCluster("Low Lambda");
    else if (clust.lambda > 1.0)
        throw InvalidCluster("High Lambda");

    // Parallax
    // Bounded between 1 and 10^5 parsecs
    else if (s.modIsParallax && clust.clust[0].mod < 0.00001)
        throw InvalidCluster("Low parallax, Clust A");

    else if (s.modIsParallax && clust.clust[1].mod < 0.00001)
        throw InvalidCluster("Low parallax, Clust B");

    else if (s.modIsParallax && clust.clust[0].mod > 1.0)
        throw InvalidCluster("High parallax, Clust A");

    else if (s.modIsParallax && clust.clust[1].mod > 1.0)
        throw InvalidCluster("High parallax, Clust B");
}

MpiMcmcApplication::MpiMcmcApplication(Settings &s,
                                       MultiPopBackingStore *mcmcStore,
                                       StarParamsBackingStore *paramsStore,
                                       FieldStarLikelihoodBackingStore *fieldStarLikelihood)
    : evoModels(makeModel(s)), settings(s)
    , gen(uint32_t(s.seed * uint32_t(2654435761))) // Applies Knuth's multiplicative hash for obfuscation (TAOCP Vol. 3)
    , mcmcStore(mcmcStore), paramsStore(paramsStore)
    , fieldStarLikelihood(fieldStarLikelihood)
    , pool(s.threads)
{
    N_WD_MASS1 = s.wdInterpolationPower == 0 ? 8000 // Old default
                                             : s.wdInterpolationPower > 0 ? 64 << s.wdInterpolationPower // 2^(5 + wdInterpolationPower)
                                                                          : 64 >> -s.wdInterpolationPower;

    ctrl.priorVar.fill(0);

    clust.lambda = settings.multiPopMcmc.lambda_start;

    clust.clust[0].feh = settings.cluster.starting.Fe_H;
    clust.clust[0].priorMean[FEH] = settings.cluster.priorMeans.Fe_H;
    ctrl.priorVar[FEH]          = settings.cluster.priorSigma.Fe_H;

    clust.clust[0].mod = settings.cluster.starting.distMod;
    clust.clust[0].priorMean[MOD] = settings.cluster.priorMeans.distMod;
    ctrl.priorVar[MOD]          = settings.cluster.priorSigma.distMod;

    clust.clust[0].abs = settings.cluster.starting.Av;
    clust.clust[0].priorMean[ABS] = fabs(settings.cluster.priorMeans.Av);
    ctrl.priorVar[ABS]          = settings.cluster.priorSigma.Av;

    clust.clust[0].age = settings.cluster.starting.logAge;
    clust.clust[0].priorMean[AGE] = settings.cluster.priorMeans.logAge;
    ctrl.priorVar[AGE]          = settings.cluster.priorSigma.logAge;

    clust.clust[0].carbonicity = settings.cluster.starting.carbonicity;
    clust.clust[0].priorMean[CARBONICITY] = settings.cluster.priorMeans.carbonicity;
    ctrl.priorVar[CARBONICITY]          = settings.cluster.priorSigma.carbonicity;

    clust.clust[0].yyy = settings.cluster.starting.Y;
    clust.clust[0].priorMean[YYY] = settings.cluster.priorMeans.Y;
    ctrl.priorVar[YYY]          = 0;

    // No IFMR code

    clust.clust[0].setM_wd_up(settings.whiteDwarf.M_wd_up);

    for (auto &var : ctrl.priorVar)
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

    std::copy(ctrl.priorVar.begin(), ctrl.priorVar.end(), clust.clust[0].priorVar.begin());

    clust.clust[1] = clust.clust[0];

    // Cluster-specific variables
    clust.clust[0].yyy = settings.multiPopMcmc.YA_start;
    clust.clust[1].yyy = settings.multiPopMcmc.YB_start;

    /* read burnIter and nIter */
    {
        ctrl.burnIter = settings.singlePopMcmc.burnIter;
        int nParamsUsed = 0;

        for (int p = 0; p < NPARAMS; p++)
        {
            if (ctrl.priorVar.at(p) > EPSILON || p == YYA || p == YYB || p == LAMBDA)
            {
                nParamsUsed++;
            }
        }

        if (ctrl.burnIter < 2 * (nParamsUsed + 1))
        {
            ctrl.burnIter = 2 * (nParamsUsed + 1);
            cerr << "burnIter below minimum allowable size. Increasing to " << ctrl.burnIter << endl;
        }
    }

    ctrl.nIter = settings.singlePopMcmc.maxIter;
    ctrl.thin = settings.singlePopMcmc.thin;
}


// Initialize SSE memory for noBinaries
// THIS LEAKS MEMORY LIKE A SIEVE DUE TO FORGETTING TALLOC
// But it's okay if we only call it once (until it gets moved to the constructor)
void MpiMcmcApplication::allocateSSEMem()
{
    using aligned_m128 = std::aligned_storage<16, 16>::type;

    if (! systems.empty())
    {
        // howManyFilts has to be a multiple of two to make the SSE code happy
        // Add 1 and round down.
        howManyFiltsAligned = ((systems.front().obsPhot.size() + 1) & ~0x1);
        howManyFilts        = systems.front().obsPhot.size();
        size_t howManyWeNeed       = systems.size() * howManyFiltsAligned;
        size_t howManyWeAlloc      = howManyWeNeed / 2;

        aligned_m128* talloc = new aligned_m128[howManyWeAlloc];
        sysVars = new(talloc) double[howManyWeNeed];

        talloc  = new aligned_m128[howManyWeAlloc];
        sysVar2 = new(talloc) double[howManyWeNeed];

        talloc = new aligned_m128[howManyWeAlloc];
        sysObs = new(talloc) double[howManyWeNeed];

    }

    int i = 0;

    for (auto s : systems)
    {
        for (size_t k = 0; k < howManyFiltsAligned; ++k, ++i)
        {
            if ((k < s.variance.size()) && (s.variance.at(k) > EPS))
            {
                sysVars[i] = s.variance.at(k) * clust.clust[0].varScale;
                sysVar2[i] = __builtin_log (TWO_M_PI * sysVars[i]);

                // Removes a division from the noBinaries loop which turns a 14 cycle loop into a 3.8 cycle loop
                sysVars[i] = 1 / sysVars[i];

                sysObs [i] = s.obsPhot.at(k);
            }
            else
            {
                sysVars[i] = 0.0;
                sysVar2[i] = 0.0;

                sysObs [i] = 0.0;
            }
        }
    }
}


int MpiMcmcApplication::run()
{
    double fsLike;

    array<double, NPARAMS> stepSize;
    std::copy(settings.singlePopMcmc.stepSize.begin(), settings.singlePopMcmc.stepSize.end(), stepSize.begin());
    stepSize.at(LAMBDA) = settings.multiPopMcmc.lambdaStep;


    if (settings.verbose)
    {
        cout << std::setprecision(3) << std::fixed;
        cout << "Model boundaries are (" << evoModels.mainSequenceEvol->getMinAge()
             << ", " << evoModels.mainSequenceEvol->getMaxAge() << ") log years." << endl;

        cout << "Binaries are " << (settings.noBinaries ? "OFF" : "ON") << endl;
    }

    if (   settings.cluster.starting.logAge < evoModels.mainSequenceEvol->getMinAge()
        || settings.cluster.starting.logAge > evoModels.mainSequenceEvol->getMaxAge())
    {
        if (! settings.overrideBounds)
        {
            cerr << std::setprecision(3) << std::fixed
                 << "\n***Error: Starting value \"logAge\" (" << settings.cluster.starting.logAge
                 << ") is outside the model boundaries (" << evoModels.mainSequenceEvol->getMinAge()
                 << ", " << evoModels.mainSequenceEvol->getMaxAge()
                 << ").\n   Use the `--overrideBounds` flag if you really want this.\n[Exiting...]"
                 << endl;

            exit(-1);
        }
        else if (settings.verbose)
        {
            cout << std::setprecision(3) << std::fixed
                 << "\n***Warning: logAge (" << settings.cluster.starting.logAge
                 << ") is outside the model boundaries (" << evoModels.mainSequenceEvol->getMinAge()
                 << ", " << evoModels.mainSequenceEvol->getMaxAge()
                 << ").\n   Continuing due to `--overrideBounds` flag." << endl;
        }
    }

    {
        // Read photometry and calculate fsLike
        vector<double> filterPriorMin;
        vector<double> filterPriorMax;

        // open files for reading (data) and writing
        // rData implcitly relies on going out of scope to close the photometry file
        // This is awful, but pretty (since this code is, at time of writing, in restricted, anonymous scope
        std::ifstream rData(settings.files.phot);

        if (!rData)
        {
            cerr << "***Error: Photometry file " << settings.files.phot << " was not found.***" << endl;
            cerr << "[Exiting...]" << endl;
            exit (-1);
        }

        auto ret = base::utility::readPhotometry (rData, filterPriorMin, filterPriorMax, settings);
        auto filterNames = ret.first;

        for (auto r : ret.second)
        {
            if (r.observedStatus == StarStatus::MSRG)
            {
                // Everything goes into the main run
                mainRunSystems.push_back(r);
                mainRunSystemNames.push_back(r.id);

                if (r.useDuringBurnIn)
                {
                    // Only burnin things go into systems
                    systems.push_back(r);
                    systemNames.push_back(r.id);
                }

                // systems must be overwritten before the main run
            }
            else if (r.observedStatus == StarStatus::WD)
                cerr << "Found unsupported WD in photometry... Continuing anyway." << endl;
            else
                cerr << "Found unsupported star in photometry, type '" << r.observedStatus << "'... Continuing anyway." << endl;
        }

        if ( systems.empty() )
        {
            cerr << "No stars loaded... Exiting." << endl;
            exit(-1);
        }

        evoModels.restrictFilters(filterNames, settings.allowInvalidModels);

        if (settings.cluster.index < 0 || static_cast<size_t>(settings.cluster.index) > filterNames.size())
        {
            cerr << "***Error: " << settings.cluster.index << " not a valid magnitude index.  Choose 0, 1,or 2.***" << endl;
            cerr << "[Exiting...]" << endl;
            exit (1);
        }

        if (systems.size() > 1)
        {
            double logFieldStarLikelihood = 0.0;

            for (size_t filt = 0; filt < filterNames.size(); filt++)
            {
                logFieldStarLikelihood -= log (filterPriorMax.at(filt) - filterPriorMin.at(filt));
            }
            fsLike = exp (logFieldStarLikelihood);
        }
        else
        {
            fsLike = 0;
        }
    }

    fieldStarLikelihood->save({fsLike});

    Chain<DualPopCluster> chain(static_cast<uint32_t>(std::uniform_int_distribution<>()(gen)), fsLike, ctrl.priorVar, clust, *mcmcStore, *paramsStore, settings.modIsParallax);

    allocateSSEMem();

    // Begin initChain
    {
        for (auto system : systems)
        {
            system.clustStarProposalDens = system.clustStarPriorDens;   // Use prior prob of being clus star
        }
    }
    // end initChain

    // Assuming fsLike doesn't change, this is the "global" logPost function
    auto logPostFunc = std::bind(&MpiMcmcApplication::logPostStep, this, _1, fsLike);

    std::function<void(const DualPopCluster&)> checkPriors = std::bind(&ensurePriors, std::cref(settings), _1);

    int adaptiveBurnIter = 0;
    bool acceptedOne = false;  // Keeps track of whether we have two accepted trials in a ros
    bool acceptedOnce = false; // Keeps track of whether we have ever accepted a trial

    // Stage 1 burnin
    {
        if ( settings.verbose )
            cout << "\nRunning Stage 1 burnin..." << flush;

        auto proposalFunc = std::bind(&MpiMcmcApplication::propClustBigSteps, this, _1, std::cref(ctrl), std::cref(stepSize));
        chain.run(AdaptiveMcmcStage::FixedBurnin, systemNames, proposalFunc, logPostFunc, checkPriors, settings.singlePopMcmc.adaptiveBigSteps);

        if ( settings.verbose )
            cout << " Complete (acceptanceRatio = " << chain.acceptanceRatio() << ")" << endl;

        chain.reset(); // Reset the chain to forget this part of the burnin.
    }

    if ( settings.verbose )
        cout << "\nRunning Stage 2 (adaptive) burnin..." << endl;

    // Run adaptive burnin (stage 2)
    // -----------------------------
    // Exits after two consecutive trials with a 20% < x < 40% acceptanceRatio,
    // or after the numeber of iterations exceeds the maximum burnin iterations.
    do
    {
        // Convenience variable
        const int trialIter = settings.singlePopMcmc.trialIter;

        // Increment the number of iterations we've gone through
        // If the step sizes aren't converging on an acceptable acceptance ratio, this can kick us out of the burnin
        adaptiveBurnIter += trialIter;

        // Reset the ratio to determine step scaling for this trial
        chain.resetRatio();

        std::function<DualPopCluster(DualPopCluster)> proposalFunc = std::bind(&MpiMcmcApplication::propClustIndep, this, _1, std::cref(ctrl), std::cref(stepSize), 5);

        if (settings.singlePopMcmc.bigStepBurnin)
        {
            // Run big steps for the entire trial
            chain.run(AdaptiveMcmcStage::AdaptiveBurnin, systemNames, proposalFunc, logPostFunc, checkPriors, trialIter);
        }
        else if (acceptedOnce)
        {
            proposalFunc = std::bind(&MpiMcmcApplication::propClustIndep, this, _1, std::cref(ctrl), std::cref(stepSize), 1);
            chain.run(AdaptiveMcmcStage::AdaptiveBurnin, systemNames, proposalFunc, logPostFunc, checkPriors, trialIter / 2);
        }
        else
        {
            // Run big steps for half the trial
            chain.run(AdaptiveMcmcStage::AdaptiveBurnin, systemNames, proposalFunc, logPostFunc, checkPriors, trialIter / 2);

            // Then run smaller steps for the second half
            proposalFunc = std::bind(&MpiMcmcApplication::propClustIndep, this, _1, std::cref(ctrl), std::cref(stepSize), 1);
            chain.run(AdaptiveMcmcStage::AdaptiveBurnin, systemNames, proposalFunc, logPostFunc, checkPriors, trialIter / 2);
        }

        double acceptanceRatio = chain.acceptanceRatio();

        if (acceptanceRatio <= 0.4 && acceptanceRatio >= 0.2)
        {
            acceptedOnce = true;

            if (acceptedOne)
            {
                if ( settings.verbose )
                    cout << "  Leaving adaptive burnin early with an acceptance ratio of " << acceptanceRatio << " (iteration " << adaptiveBurnIter + settings.singlePopMcmc.adaptiveBigSteps << ")" << endl;

                break;
            }
            else
            {
                if ( settings.verbose )
                    cout << "    Acceptance ratio: " << acceptanceRatio << ". Trying for trend." << endl;
                acceptedOne = true;
            }
        }
        else
        {
            acceptedOne = false;

            if ( settings.verbose )
                cout << "    Acceptance ratio: " << acceptanceRatio << ". Retrying." << endl;

            scaleStepSizes(stepSize, chain.acceptanceRatio()); // Adjust step sizes
        }
    } while (adaptiveBurnIter < ctrl.burnIter);


    // Main run
    // Overwrite the burnin photometry set with the entire photometry set
    systems = mainRunSystems;
    systemNames = mainRunSystemNames;

    // Don't forget to repopulate the SSE memory
    allocateSSEMem();

    {
        // Make sure and pull the covariance matrix before resetting the chain
        auto proposalFunc = std::bind(&MpiMcmcApplication::propClustCorrelated, this, _1, std::cref(ctrl), chain.makeCholeskyDecomp());

        chain.reset();

        if ( settings.verbose )
            cout << "\nStarting adaptive run... " << flush;

        chain.run(AdaptiveMcmcStage::AdaptiveMainRun, systemNames, proposalFunc, logPostFunc, checkPriors, settings.singlePopMcmc.stage3Iter);

        if ( settings.verbose )
            cout << " Preliminary acceptanceRatio = " << chain.acceptanceRatio() << endl;
    }

    // Begin main run
    // Main run proceeds in increments of 1, adapting the covariance matrix after every increment
    for (auto iters = 0; iters < ctrl.nIter; ++iters)
    {
        auto proposalFunc = std::bind(&MpiMcmcApplication::propClustCorrelated, this, _1, std::cref(ctrl), chain.makeCholeskyDecomp());

        chain.run(AdaptiveMcmcStage::MainRun, systemNames, proposalFunc, logPostFunc, checkPriors, 1, ctrl.thin);
    }

    if ( settings.verbose )
    {
        cout << "\nFinal acceptance ratio: " << chain.acceptanceRatio() << "\n" << endl;

        MsBoundsError::countSummary();
    }

    return 0;
}


void MpiMcmcApplication::scaleStepSizes (array<double, NPARAMS> &stepSize, double acceptanceRatio)
{
    function<double(double)> scaleFactor = [](double acceptanceRatio) {
        double factor = 1.0;

        if (acceptanceRatio > 0.9)
        {
            factor = 2.0;
        }
        else if (acceptanceRatio > 0.7)
        {
            factor = 1.8;
        }
        else if (acceptanceRatio > 0.5)
        {
            factor = 1.5;
        }
        else if (acceptanceRatio > 0.4)
        {
            factor = 1.2;
        }
        else if (acceptanceRatio < 0.2)
        {
            factor = 1 / 1.5;
        }
        else if (acceptanceRatio < 0.15)
        {
            factor = 1 / 1.8;
        }
        else if (acceptanceRatio < 0.05)
        {
            factor = 0.5;
        }

        return factor;
    };

    for (int p = 0; p < NPARAMS; p++)
    {
        if (ctrl.priorVar.at(p) > EPSILON || p == YYA || p == YYB || p == LAMBDA)
        {
            stepSize.at(p) *= scaleFactor(acceptanceRatio);
        }
    }
}


DualPopCluster MpiMcmcApplication::propClustBigSteps (DualPopCluster clust, struct ifmrMcmcControl const &ctrl, array<double, NPARAMS> const &stepSize)
{
    return propClustIndep(clust, ctrl, stepSize, 25.0);
}

DualPopCluster MpiMcmcApplication::propClustIndep (DualPopCluster clust, struct ifmrMcmcControl const &ctrl, array<double, NPARAMS> const &stepSize, double scale)
{
    /* DOF defined in densities.h */
    int p;

    for (p = 0; p < NPARAMS; p++)
    {
        if (ctrl.priorVar.at(p) > EPSILON && p != YYA && p != YYB && p != LAMBDA)
        {
            clust.clust[0].setParam(p, clust.clust[0].getParam(p) + sampleT (gen, scale * stepSize.at(p) * stepSize.at(p)));
            clust.clust[1].setParam(p, clust.clust[0].getParam(p));
        }
    }

    clust.clust[0].yyy += sampleT (gen, scale * stepSize.at(YYY) * stepSize.at(YYY));
    clust.clust[1].yyy += sampleT (gen, scale * stepSize.at(YYY) * stepSize.at(YYY));

    clust.lambda += sampleT (gen, scale * stepSize.at(LAMBDA) * stepSize.at(LAMBDA));

    return clust;
}

DualPopCluster MpiMcmcApplication::propClustCorrelated (DualPopCluster clust, struct ifmrMcmcControl const &ctrl, Matrix<double, NPARAMS, NPARAMS> const &propMatrix)
{
    array<double, NPARAMS> tDraws;

    for (auto &d : tDraws)
    {
        d = sampleT (gen, 1.0);
    }

    for (int p = 0; p < NPARAMS; ++p)
    {
        if (ctrl.priorVar.at(p) > EPSILON || p == YYA || p == YYB || p == LAMBDA)
        {
            double corrProps = 0;

            for (int k = 0; k <= p; ++k)
            {                           /* propMatrix is lower diagonal */
                if (ctrl.priorVar.at(k) > EPSILON || k == YYA || k == YYB || k == LAMBDA)
                {
                    corrProps += propMatrix.at(p).at(k) * tDraws[k];
                }
            }

            if (p == YYA)
                clust.clust[0].yyy += corrProps;
            else if (p == YYB)
                clust.clust[1].yyy += corrProps;
            else if (p == LAMBDA)
                clust.lambda += corrProps;
            else
            {
                // Set the value in clust[0]
                clust.clust[0].setParam(p, clust.clust[0].getParam(p) + corrProps);
                // Then copy the value from clust[0] to clust[1]
                clust.clust[1].setParam(p, clust.clust[0].getParam(p));
            }
        }
    }

    return clust;
}

tuple<double, vector<double>> MpiMcmcApplication::logPostStep(DualPopCluster &propClust, double fsLike)
{
    double logPostProp = propClust.clust[0].logPrior (evoModels);

    array<unique_ptr<Isochrone>, 2> isochrone;
    array<Cluster*, 2> cluster = { &propClust.clust[0], &propClust.clust[1] };

    vector<double> starData;

    pool.parallelFor(2, [=,&isochrone](int i)
    {
        isochrone[i].reset(evoModels.mainSequenceEvol->deriveIsochrone(cluster[i]->feh, cluster[i]->yyy, cluster[i]->age));
    });

    auto sSize = systems.size();

    const double lambdaA = propClust.lambda;
    const double lambdaB = (1 - lambdaA);

    /* marginalize over isochrone */
    array<vector<double>, 2> post;
    auto &postA = post[0];
    auto &postB = post[1];

    pool.parallelFor(2, [=, &isochrone, &propClust, &post](int i)
    {
        if (settings.noBinaries)
        {
            post[i] = margEvolveNoBinaries (propClust.clust[i], evoModels, *isochrone.at(i), pool, sysVars, sysVar2, sysObs, sSize, howManyFiltsAligned, howManyFilts, settings.modIsParallax, settings.eepInterpolationPower);
        }
        else
        {
            post[i] = margEvolveWithBinary (propClust.clust[i], systems, evoModels, *isochrone.at(i), pool, settings.modIsParallax, settings.eepInterpolationPower);
        }
    });


    for (size_t i = 0; i < sSize; ++i)
    {
        starData.push_back(postA[i]);
        starData.push_back(postB[i]);

        // Give both A and B a chance to be the only contributor for this star.
        if (postA[i] > 0.0)
        {
            postA[i] *= lambdaA;
        }
        else
        {
            postA[i] = 0.0;
        }

        // Independent of above
        if (postB[i] > 0.0)
        {
            postA[i] += lambdaB * postB[i];
        }

        // It's also possible that it comes out zero and we waste our time doing this multiplication
        // But it's faster to do it than to check for zero.
        postA[i] *= systems[i].clustStarPriorDens;

        /* marginalize over field star status */
        logPostProp += __builtin_log ((1.0 - systems[i].clustStarPriorDens) * fsLike + postA[i]);
    }

    return std::make_tuple(logPostProp, starData);
}


double MpiMcmcApplication::wdGridMass (int point) const
{
    // I think this only gets calculated once, but I can't quite be sure...
    // It may not even matter, but every little bit helps, eh? And hey, preoptimization...
    static const double dWdMass1 = (settings.whiteDwarf.M_wd_up - MIN_MASS1) / (double) N_WD_MASS1;

    return MIN_MASS1 + (dWdMass1 * point);
}
