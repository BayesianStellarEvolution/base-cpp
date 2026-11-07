// Copyright (c) Embry-Riddle Aeronautical University. All rights reserved.

#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <memory>
#include <random>
#include <string>
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
    static std::string photFile;

protected:
    static Model* evoModels;
    static Isochrone* isochrone;
    static Cluster clust;
    static size_t howManyFilts;
    static size_t howManyFiltsAligned;
    static bool initialized;
    static std::vector<StellarSystem> loadedSystems;

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

        // If photFile provided, load and use its filters; otherwise use test filters
        if (!photFile.empty()) {
            std::ifstream fin(photFile);
            if (!fin) {
                throw std::runtime_error("Cannot open photometry file: " + photFile);
            }
            std::vector<double> filterPriorMin, filterPriorMax;
            auto [filterNames, systems] = base::utility::readPhotometry(fin, filterPriorMin, filterPriorMax, settings);
            evoModels->restrictFilters(filterNames, false);
            loadedSystems = std::move(systems);
            std::cerr << "Loaded " << loadedSystems.size() << " stars from " << photFile << "\n";
        } else {
            std::vector<std::string> testFilters = {"G", "G_BP", "G_RP"};
            evoModels->restrictFilters(testFilters, false);
        }

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

    // Returns loaded systems if photFile was provided, otherwise generates n synthetic systems
    // If n < loadedSystems.size(), returns first n loaded systems
    vector<StellarSystem> getSystems(size_t n = 100, unsigned seed = 42) {
        if (!loadedSystems.empty()) {
            if (n >= loadedSystems.size()) {
                return loadedSystems;
            }
            return vector<StellarSystem>(loadedSystems.begin(), loadedSystems.begin() + n);
        }
        return createSystems(n, seed);
    }

    // Prepare aligned arrays for SIMD (both x86 SSE and ARM NEON)
    void prepareAlignedArrays(const vector<StellarSystem>& systems,
                              vector<double>& sysVars,
                              vector<double>& sysVar2,
                              vector<double>& sysObs,
                              vector<double>& sysFluxLo,
                              vector<double>& sysFluxHi) {
        size_t n = systems.size();
        size_t totalSize = n * howManyFiltsAligned;

        sysVars.resize(totalSize);
        sysVar2.resize(totalSize);
        sysObs.resize(totalSize);
        sysFluxLo.resize(totalSize);
        sysFluxHi.resize(totalSize);

        const double TWO_M_PI = 2 * M_PI;
        const double nSD = 3.0;
        size_t idx = 0;
        for (const auto& sys : systems) {
            for (size_t k = 0; k < howManyFiltsAligned; ++k, ++idx) {
                if (k < howManyFilts) {
                    double var = sys.variance[k] * clust.varScale;
                    sysVars[idx] = 1.0 / var;
                    sysVar2[idx] = std::log(TWO_M_PI * var);
                    sysObs[idx] = sys.obsPhot[k];
                    double nSDTerm = nSD * std::sqrt(var);
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
std::string NoBinariesTest::photFile;
std::vector<StellarSystem> NoBinariesTest::loadedSystems;

TEST_F(NoBinariesTest, ReturnsFiniteValues) {
    ASSERT_TRUE(initialized) << "Models not initialized";

    auto systems = getSystems(100);
    vector<double> sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi;
    prepareAlignedArrays(systems, sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi);

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

    auto systems = getSystems(50, 123);
    vector<double> sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi;
    prepareAlignedArrays(systems, sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi);

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

    auto systems = getSystems(50, 456);

    // Get reference results (always uses scalar implementation)
    auto refResults = margEvolveNoBinaries_Reference(
        clust, *evoModels, *isochrone,
        systems, false, 0);

    // Compare SIMD (x86 SSE or ARM NEON) against scalar reference
    vector<double> sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi;
    prepareAlignedArrays(systems, sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi);

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

// Test NoBinaries with single system
TEST_F(NoBinariesTest, NoBinariesSingleSystem) {
    ASSERT_TRUE(initialized) << "Models not initialized";

    auto systems = getSystems(1, 777);

    auto refResults = margEvolveNoBinaries_Reference(
        clust, *evoModels, *isochrone,
        systems, false, 0);

    vector<double> sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi;
    prepareAlignedArrays(systems, sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi);

    auto simdResults = margEvolveNoBinaries(
        clust, *evoModels, *isochrone, pool,
        sysVars.data(), sysVar2.data(), sysObs.data(),
        systems.size(), howManyFiltsAligned, howManyFilts,
        false, 0);

    ASSERT_EQ(simdResults.size(), 1u);
    ASSERT_EQ(refResults.size(), 1u);

    if (refResults[0] != 0.0) {
        double relError = std::abs((simdResults[0] - refResults[0]) / refResults[0]);
        EXPECT_LT(relError, 1e-10);
    }
}

// Test NoBinaries with larger count
TEST_F(NoBinariesTest, NoBinariesLargeCount) {
    ASSERT_TRUE(initialized) << "Models not initialized";

    auto systems = getSystems(500, 888);

    auto refResults = margEvolveNoBinaries_Reference(
        clust, *evoModels, *isochrone,
        systems, false, 0);

    vector<double> sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi;
    prepareAlignedArrays(systems, sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi);

    auto simdResults = margEvolveNoBinaries(
        clust, *evoModels, *isochrone, pool,
        sysVars.data(), sysVar2.data(), sysObs.data(),
        systems.size(), howManyFiltsAligned, howManyFilts,
        false, 0);

    ASSERT_EQ(simdResults.size(), refResults.size());

    double maxRelError = 0.0;
    size_t mismatchCount = 0;
    for (size_t i = 0; i < simdResults.size(); ++i) {
        if (refResults[i] != 0.0) {
            double relError = std::abs((simdResults[i] - refResults[i]) / refResults[i]);
            maxRelError = std::max(maxRelError, relError);
            if (relError >= 1e-10) mismatchCount++;
        }
    }
    EXPECT_EQ(mismatchCount, 0u) << "Found " << mismatchCount << " mismatches in 500 systems";
    std::cerr << "NoBinaries large count max relative error: " << maxRelError << "\n";
}

TEST_F(NoBinariesTest, BinarySIMDMatchesReference) {
    ASSERT_TRUE(initialized) << "Models not initialized";

    auto systems = getSystems(50, 999);

    // Get reference results
    auto refResults = margEvolveWithBinary(
        clust, systems, *evoModels, *isochrone, pool,
        false, 0);

    // Prepare aligned arrays
    vector<double> sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi;
    prepareAlignedArrays(systems, sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi);

    // Get SIMD results
    auto simdResults = margEvolveWithBinary(
        clust, *evoModels, *isochrone, pool,
        systems,
        sysVars.data(), sysVar2.data(), sysObs.data(),
        sysFluxLo.data(), sysFluxHi.data(),
        howManyFiltsAligned, howManyFilts,
        false, 0);

    ASSERT_EQ(simdResults.size(), refResults.size());

    double maxRelError = 0.0;
    for (size_t i = 0; i < simdResults.size(); ++i) {
        if (refResults[i] != 0.0) {
            double relError = std::abs((simdResults[i] - refResults[i]) / refResults[i]);
            maxRelError = std::max(maxRelError, relError);
            EXPECT_LT(relError, 1e-10)
                << "Result[" << i << "]: simd=" << simdResults[i]
                << " reference=" << refResults[i]
                << " relError=" << relError;
        }
    }
    std::cerr << "Binary SIMD max relative error: " << maxRelError << "\n";
}

// Test with single system (edge case)
TEST_F(NoBinariesTest, BinarySIMDSingleSystem) {
    ASSERT_TRUE(initialized) << "Models not initialized";

    auto systems = getSystems(1, 111);

    auto refResults = margEvolveWithBinary(
        clust, systems, *evoModels, *isochrone, pool,
        false, 0);

    vector<double> sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi;
    prepareAlignedArrays(systems, sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi);

    auto simdResults = margEvolveWithBinary(
        clust, *evoModels, *isochrone, pool,
        systems,
        sysVars.data(), sysVar2.data(), sysObs.data(),
        sysFluxLo.data(), sysFluxHi.data(),
        howManyFiltsAligned, howManyFilts,
        false, 0);

    ASSERT_EQ(simdResults.size(), 1u);
    ASSERT_EQ(refResults.size(), 1u);

    if (refResults[0] != 0.0) {
        double relError = std::abs((simdResults[0] - refResults[0]) / refResults[0]);
        EXPECT_LT(relError, 1e-10);
    }
}

// Test determinism for binary SIMD
TEST_F(NoBinariesTest, BinarySIMDDeterministic) {
    ASSERT_TRUE(initialized) << "Models not initialized";

    auto systems = getSystems(25, 222);
    vector<double> sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi;
    prepareAlignedArrays(systems, sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi);

    auto results1 = margEvolveWithBinary(
        clust, *evoModels, *isochrone, pool,
        systems,
        sysVars.data(), sysVar2.data(), sysObs.data(),
        sysFluxLo.data(), sysFluxHi.data(),
        howManyFiltsAligned, howManyFilts,
        false, 0);

    auto results2 = margEvolveWithBinary(
        clust, *evoModels, *isochrone, pool,
        systems,
        sysVars.data(), sysVar2.data(), sysObs.data(),
        sysFluxLo.data(), sysFluxHi.data(),
        howManyFiltsAligned, howManyFilts,
        false, 0);

    ASSERT_EQ(results1.size(), results2.size());
    for (size_t i = 0; i < results1.size(); ++i) {
        EXPECT_DOUBLE_EQ(results1[i], results2[i])
            << "Binary SIMD results differ at index " << i;
    }
}

// Test with larger system count
TEST_F(NoBinariesTest, BinarySIMDLargeCount) {
    ASSERT_TRUE(initialized) << "Models not initialized";

    auto systems = getSystems(500, 333);

    auto refResults = margEvolveWithBinary(
        clust, systems, *evoModels, *isochrone, pool,
        false, 0);

    vector<double> sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi;
    prepareAlignedArrays(systems, sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi);

    auto simdResults = margEvolveWithBinary(
        clust, *evoModels, *isochrone, pool,
        systems,
        sysVars.data(), sysVar2.data(), sysObs.data(),
        sysFluxLo.data(), sysFluxHi.data(),
        howManyFiltsAligned, howManyFilts,
        false, 0);

    ASSERT_EQ(simdResults.size(), refResults.size());

    double maxRelError = 0.0;
    size_t mismatchCount = 0;
    for (size_t i = 0; i < simdResults.size(); ++i) {
        if (refResults[i] != 0.0) {
            double relError = std::abs((simdResults[i] - refResults[i]) / refResults[i]);
            maxRelError = std::max(maxRelError, relError);
            if (relError >= 1e-10) mismatchCount++;
        }
    }
    EXPECT_EQ(mismatchCount, 0u) << "Found " << mismatchCount << " mismatches in 500 systems";
    std::cerr << "Large count max relative error: " << maxRelError << "\n";
}

// Test with extreme mass values near isochrone boundaries
TEST_F(NoBinariesTest, BinarySIMDExtremeMasses) {
    ASSERT_TRUE(initialized) << "Models not initialized";

    // Create systems with masses at extreme ends of the isochrone
    vector<StellarSystem> systems;
    std::mt19937 gen(444);
    std::uniform_real_distribution<> varDist(0.001, 0.1);
    std::uniform_real_distribution<> obsDist(10.0, 20.0);

    // Get mass range from isochrone
    double minMass = isochrone->eeps.front().mass;
    double maxMass = isochrone->eeps.back().mass;

    // Test masses near boundaries
    std::vector<double> testMasses = {
        minMass + 0.001,          // Just above minimum
        minMass + 0.1,            // Near minimum
        (minMass + maxMass) / 2,  // Middle
        maxMass - 0.1,            // Near maximum
        maxMass - 0.001           // Just below maximum
    };

    for (double mass : testMasses) {
        StellarSystem sys;
        sys.primary.mass = mass;
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

    auto refResults = margEvolveWithBinary(
        clust, systems, *evoModels, *isochrone, pool,
        false, 0);

    vector<double> sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi;
    prepareAlignedArrays(systems, sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi);

    auto simdResults = margEvolveWithBinary(
        clust, *evoModels, *isochrone, pool,
        systems,
        sysVars.data(), sysVar2.data(), sysObs.data(),
        sysFluxLo.data(), sysFluxHi.data(),
        howManyFiltsAligned, howManyFilts,
        false, 0);

    ASSERT_EQ(simdResults.size(), refResults.size());

    for (size_t i = 0; i < simdResults.size(); ++i) {
        EXPECT_TRUE(std::isfinite(simdResults[i]))
            << "SIMD result[" << i << "] (mass=" << testMasses[i] << ") is not finite";
        EXPECT_TRUE(std::isfinite(refResults[i]))
            << "Reference result[" << i << "] (mass=" << testMasses[i] << ") is not finite";

        if (refResults[i] != 0.0) {
            double relError = std::abs((simdResults[i] - refResults[i]) / refResults[i]);
            EXPECT_LT(relError, 1e-10)
                << "Mismatch at mass=" << testMasses[i];
        }
    }
}

// Test with varying variances (small and large)
TEST_F(NoBinariesTest, BinarySIMDVarianceRange) {
    ASSERT_TRUE(initialized) << "Models not initialized";

    vector<StellarSystem> systems;
    std::mt19937 gen(555);
    std::uniform_real_distribution<> massDist(0.5, 1.2);
    std::uniform_real_distribution<> obsDist(12.0, 18.0);

    // Test different variance scales
    std::vector<double> varScales = {0.0001, 0.001, 0.01, 0.1, 1.0};

    for (double varScale : varScales) {
        StellarSystem sys;
        sys.primary.mass = massDist(gen);
        sys.secondary.mass = 0.0;
        sys.clustStarPriorDens = 1.0;
        sys.obsPhot.resize(howManyFilts);
        sys.variance.resize(howManyFilts);
        for (size_t f = 0; f < howManyFilts; ++f) {
            sys.obsPhot[f] = obsDist(gen);
            sys.variance[f] = varScale;
        }
        systems.push_back(sys);
    }

    auto refResults = margEvolveWithBinary(
        clust, systems, *evoModels, *isochrone, pool,
        false, 0);

    vector<double> sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi;
    prepareAlignedArrays(systems, sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi);

    auto simdResults = margEvolveWithBinary(
        clust, *evoModels, *isochrone, pool,
        systems,
        sysVars.data(), sysVar2.data(), sysObs.data(),
        sysFluxLo.data(), sysFluxHi.data(),
        howManyFiltsAligned, howManyFilts,
        false, 0);

    ASSERT_EQ(simdResults.size(), refResults.size());

    for (size_t i = 0; i < simdResults.size(); ++i) {
        EXPECT_TRUE(std::isfinite(simdResults[i]))
            << "SIMD result[" << i << "] (var=" << varScales[i] << ") is not finite";

        if (refResults[i] != 0.0 && std::isfinite(refResults[i])) {
            double relError = std::abs((simdResults[i] - refResults[i]) / refResults[i]);
            EXPECT_LT(relError, 1e-10)
                << "Mismatch at variance scale=" << varScales[i];
        }
    }
}

// Test with different EEP interpolation powers
TEST_F(NoBinariesTest, BinarySIMDInterpolationPowers) {
    ASSERT_TRUE(initialized) << "Models not initialized";

    auto systems = getSystems(20, 666);
    vector<double> sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi;
    prepareAlignedArrays(systems, sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi);

    // Test different interpolation powers (-3 to 3)
    for (int power = -3; power <= 3; ++power) {
        auto refResults = margEvolveWithBinary(
            clust, systems, *evoModels, *isochrone, pool,
            false, power);

        auto simdResults = margEvolveWithBinary(
            clust, *evoModels, *isochrone, pool,
            systems,
            sysVars.data(), sysVar2.data(), sysObs.data(),
            sysFluxLo.data(), sysFluxHi.data(),
            howManyFiltsAligned, howManyFilts,
            false, power);

        ASSERT_EQ(simdResults.size(), refResults.size());

        double maxRelError = 0.0;
        for (size_t i = 0; i < simdResults.size(); ++i) {
            if (refResults[i] != 0.0) {
                double relError = std::abs((simdResults[i] - refResults[i]) / refResults[i]);
                maxRelError = std::max(maxRelError, relError);
                EXPECT_LT(relError, 1e-10)
                    << "Mismatch at power=" << power << ", index=" << i;
            }
        }
        std::cerr << "EEP power " << power << " max relative error: " << maxRelError << "\n";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    // Parse --photFile before passing to test suite
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--photFile" && i + 1 < argc) {
            NoBinariesTest::photFile = argv[i + 1];
            break;
        }
    }

    // Pass remaining args to the test suite for config loading
    NoBinariesTest::savedArgc = argc;
    NoBinariesTest::savedArgv = argv;

    return RUN_ALL_TESTS();
}
