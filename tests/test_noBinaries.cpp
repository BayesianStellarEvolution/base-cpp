// Copyright (c) Embry-Riddle Aeronautical University. All rights reserved.

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <random>
#include <vector>

#include "Cluster.hpp"
#include "Filters.hpp"
#include "GlobalSettings.hpp"
#include "Model.hpp"
#include "Star.hpp"
#include "Utility.hpp"
#include "marg.hpp"

using std::vector;

class NoBinariesTest : public ::testing::Test {
public:
    static int savedArgc;
    static char** savedArgv;

protected:
    static Model* evoModels;
    static Isochrone* isochrone;
    static Cluster clust;
    static size_t howManyFilts;
    static size_t howManyFiltsAligned;
    static bool initialized;

    static void SetUpTestSuite() {
        if (initialized) return;

        GlobalSettings settings;
        settings.loadSettings(savedArgc, savedArgv);

        clust.feh = settings.cluster.starting.Fe_H;
        clust.mod = settings.cluster.starting.distMod;
        clust.abs = settings.cluster.starting.Av;
        clust.age = settings.cluster.starting.logAge;
        clust.carbonicity = settings.cluster.starting.carbonicity;
        clust.yyy = settings.cluster.starting.Y;
        clust.ifmrSlope = 0.08;
        clust.ifmrIntercept = 0.65;
        clust.ifmrQuadCoef = 0.0001;
        clust.setM_wd_up(settings.whiteDwarf.M_wd_up);
        clust.varScale = 1.0;

        clust.priorMean[FEH] = settings.cluster.priorMeans.Fe_H;
        clust.priorMean[MOD] = settings.cluster.priorMeans.distMod;
        clust.priorMean[ABS] = settings.cluster.priorMeans.Av;
        clust.priorMean[AGE] = settings.cluster.priorMeans.logAge;
        clust.priorMean[CARBONICITY] = settings.cluster.priorMeans.carbonicity;
        clust.priorMean[YYY] = settings.cluster.priorMeans.Y;
        clust.priorVar.fill(0.01);

        evoModels = new Model(makeModel(settings));

        // Restrict to test filters (Gaia DR3)
        std::vector<std::string> testFilters = {"G", "G_BP", "G_RP"};
        evoModels->restrictFilters(testFilters, false);

        isochrone = evoModels->mainSequenceEvol->deriveIsochrone(clust.feh, clust.yyy, clust.age);
        howManyFilts = evoModels->absCoeffs.size();
        howManyFiltsAligned = ((howManyFilts + 1) & ~0x1);

        initialized = true;
    }

    static void TearDownTestSuite() {
        delete evoModels;
        delete isochrone;
        evoModels = nullptr;
        isochrone = nullptr;
    }

    vector<StellarSystem> createSystems(size_t n, unsigned seed = 42) {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<> massDist(0.2, 1.5);
        std::uniform_real_distribution<> varDist(0.001, 0.1);
        std::uniform_real_distribution<> obsDist(10.0, 20.0);

        vector<StellarSystem> systems;
        systems.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            StellarSystem sys;
            sys.primary.mass = massDist(gen);
            sys.secondary.mass = 0.0;
            sys.clustStarPriorDens = 1.0;

            sys.obsPhot.resize(howManyFilts);
            sys.variance.resize(howManyFilts);

            for (size_t f = 0; f < howManyFilts; ++f) {
                sys.obsPhot[f] = obsDist(gen);
                sys.variance[f] = varDist(gen);
            }

            systems.push_back(sys);
        }
        return systems;
    }

    // Prepare aligned arrays for SIMD (both x86 SSE and ARM NEON)
    void prepareAlignedArrays(const vector<StellarSystem>& systems,
                              vector<double>& sysVars,
                              vector<double>& sysVar2,
                              vector<double>& sysObs) {
        size_t n = systems.size();
        size_t totalSize = n * howManyFiltsAligned;

        sysVars.resize(totalSize);
        sysVar2.resize(totalSize);
        sysObs.resize(totalSize);

        const double TWO_M_PI = 2 * M_PI;
        size_t idx = 0;
        for (const auto& sys : systems) {
            for (size_t k = 0; k < howManyFiltsAligned; ++k, ++idx) {
                if (k < howManyFilts) {
                    double var = sys.variance[k] * clust.varScale;
                    sysVars[idx] = 1.0 / var;
                    sysVar2[idx] = std::log(TWO_M_PI * var);
                    sysObs[idx] = sys.obsPhot[k];
                } else {
                    sysVars[idx] = 0.0;
                    sysVar2[idx] = 0.0;
                    sysObs[idx] = 0.0;
                }
            }
        }
    }

    base::utility::ThreadPool pool{1};
};

Model* NoBinariesTest::evoModels = nullptr;
Isochrone* NoBinariesTest::isochrone = nullptr;
Cluster NoBinariesTest::clust;
size_t NoBinariesTest::howManyFilts = 0;
size_t NoBinariesTest::howManyFiltsAligned = 0;
bool NoBinariesTest::initialized = false;
int NoBinariesTest::savedArgc = 0;
char** NoBinariesTest::savedArgv = nullptr;

TEST_F(NoBinariesTest, ReturnsFiniteValues) {
    ASSERT_TRUE(initialized) << "Models not initialized";

    auto systems = createSystems(100);
    vector<double> sysVars, sysVar2, sysObs;
    prepareAlignedArrays(systems, sysVars, sysVar2, sysObs);

    auto results = margEvolveNoBinaries(
        clust, *evoModels, *isochrone, pool,
        sysVars.data(), sysVar2.data(), sysObs.data(),
        systems.size(), howManyFiltsAligned, howManyFilts,
        false, 0);

    ASSERT_EQ(results.size(), systems.size());
    for (size_t i = 0; i < results.size(); ++i) {
        EXPECT_TRUE(std::isfinite(results[i]))
            << "Result[" << i << "] is not finite: " << results[i];
    }
}

TEST_F(NoBinariesTest, DeterministicResults) {
    ASSERT_TRUE(initialized) << "Models not initialized";

    auto systems = createSystems(50, 123);
    vector<double> sysVars, sysVar2, sysObs;
    prepareAlignedArrays(systems, sysVars, sysVar2, sysObs);

    auto results1 = margEvolveNoBinaries(
        clust, *evoModels, *isochrone, pool,
        sysVars.data(), sysVar2.data(), sysObs.data(),
        systems.size(), howManyFiltsAligned, howManyFilts,
        false, 0);

    auto results2 = margEvolveNoBinaries(
        clust, *evoModels, *isochrone, pool,
        sysVars.data(), sysVar2.data(), sysObs.data(),
        systems.size(), howManyFiltsAligned, howManyFilts,
        false, 0);

    ASSERT_EQ(results1.size(), results2.size());
    for (size_t i = 0; i < results1.size(); ++i) {
        EXPECT_DOUBLE_EQ(results1[i], results2[i])
            << "Results differ at index " << i;
    }
}

TEST_F(NoBinariesTest, HandlesEmptyInput) {
    ASSERT_TRUE(initialized) << "Models not initialized";

    vector<double> sysVars, sysVar2, sysObs;
    auto results = margEvolveNoBinaries(
        clust, *evoModels, *isochrone, pool,
        sysVars.data(), sysVar2.data(), sysObs.data(),
        0, howManyFiltsAligned, howManyFilts,
        false, 0);

    EXPECT_TRUE(results.empty());
}

TEST_F(NoBinariesTest, MatchesScalarReference) {
    ASSERT_TRUE(initialized) << "Models not initialized";

    auto systems = createSystems(50, 456);

    // Get reference results (always uses scalar implementation)
    auto refResults = margEvolveNoBinaries_Reference(
        clust, *evoModels, *isochrone,
        systems, false, 0);

    // Compare SIMD (x86 SSE or ARM NEON) against scalar reference
    vector<double> sysVars, sysVar2, sysObs;
    prepareAlignedArrays(systems, sysVars, sysVar2, sysObs);

    auto optimizedResults = margEvolveNoBinaries(
        clust, *evoModels, *isochrone, pool,
        sysVars.data(), sysVar2.data(), sysObs.data(),
        systems.size(), howManyFiltsAligned, howManyFilts,
        false, 0);

    ASSERT_EQ(optimizedResults.size(), refResults.size());

    double maxRelError = 0.0;
    for (size_t i = 0; i < optimizedResults.size(); ++i) {
        if (refResults[i] != 0.0) {
            double relError = std::abs((optimizedResults[i] - refResults[i]) / refResults[i]);
            maxRelError = std::max(maxRelError, relError);

            // Allow small relative error due to SLEEF vs standard exp
            EXPECT_LT(relError, 1e-10)
                << "Result[" << i << "]: optimized=" << optimizedResults[i]
                << " reference=" << refResults[i]
                << " relError=" << relError;
        } else {
            EXPECT_DOUBLE_EQ(optimizedResults[i], refResults[i])
                << "Result[" << i << "] differs when reference is zero";
        }
    }

    std::cerr << "Max relative error: " << maxRelError << "\n";
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    // Pass remaining args to the test suite for config loading
    NoBinariesTest::savedArgc = argc;
    NoBinariesTest::savedArgv = argv;

    return RUN_ALL_TESTS();
}
