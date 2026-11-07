#include <benchmark/benchmark.h>

#include <memory>
#include <random>
#include <vector>
#include <iostream>

#include "Cluster.hpp"
#include "Filters.hpp"
#include "GlobalSettings.hpp"
#include "Model.hpp"
#include "Star.hpp"
#include "Utility.hpp"
#include "marg.hpp"

using std::vector;

// Global state initialized once from config
static Model* g_evoModels = nullptr;
static Isochrone* g_isochrone = nullptr;
static Cluster g_clust;
static size_t g_howManyFilts = 0;
static size_t g_howManyFiltsAligned = 0;
static bool g_initialized = false;

void initializeFromSettings(int argc, char** argv) {
    if (g_initialized) return;

    GlobalSettings settings;
    settings.loadSettings(argc, argv);

    // Set up cluster parameters from config
    g_clust.feh = settings.cluster.starting.Fe_H;
    g_clust.mod = settings.cluster.starting.distMod;
    g_clust.abs = settings.cluster.starting.Av;
    g_clust.age = settings.cluster.starting.logAge;
    g_clust.carbonicity = settings.cluster.starting.carbonicity;
    g_clust.yyy = settings.cluster.starting.Y;
    g_clust.ifmrSlope = 0.08;
    g_clust.ifmrIntercept = 0.65;
    g_clust.ifmrQuadCoef = 0.0001;
    g_clust.setM_wd_up(settings.whiteDwarf.M_wd_up);
    g_clust.varScale = 1.0;

    g_clust.priorMean[FEH] = settings.cluster.priorMeans.Fe_H;
    g_clust.priorMean[MOD] = settings.cluster.priorMeans.distMod;
    g_clust.priorMean[ABS] = settings.cluster.priorMeans.Av;
    g_clust.priorMean[AGE] = settings.cluster.priorMeans.logAge;
    g_clust.priorMean[CARBONICITY] = settings.cluster.priorMeans.carbonicity;
    g_clust.priorMean[YYY] = settings.cluster.priorMeans.Y;
    g_clust.priorVar.fill(0.01);

    // Create models from settings
    g_evoModels = new Model(makeModel(settings));

    // Restrict to a typical filter set for benchmarking (Gaia DR3)
    std::vector<std::string> benchFilters = {"G", "G_BP", "G_RP"};
    g_evoModels->restrictFilters(benchFilters, false);

    // Generate isochrone
    g_isochrone = g_evoModels->mainSequenceEvol->deriveIsochrone(g_clust.feh, g_clust.yyy, g_clust.age);

    // Get filter count from model
    g_howManyFilts = g_evoModels->absCoeffs.size();
    g_howManyFiltsAligned = ((g_howManyFilts + 1) & ~0x1);

    g_initialized = true;

    std::cerr << "Initialized with " << g_howManyFilts << " filters, "
              << g_isochrone->eeps.size() << " EEPs\n";
}

class NoBinariesFixture : public benchmark::Fixture {
public:
    void SetUp(const benchmark::State& state) override {
        if (!g_initialized) {
            throw std::runtime_error("Models not initialized");
        }

        nSystems = state.range(0);

        std::mt19937 gen(42);
        std::uniform_real_distribution<> massDist(0.2, 1.5);
        std::uniform_real_distribution<> varDist(0.001, 0.1);
        std::uniform_real_distribution<> obsDist(10.0, 20.0);

        systems.clear();
        systems.reserve(nSystems);
        for (size_t i = 0; i < nSystems; ++i) {
            StellarSystem sys;
            sys.primary.mass = massDist(gen);
            sys.secondary.mass = 0.0;
            sys.clustStarPriorDens = 1.0;

            sys.obsPhot.resize(g_howManyFilts);
            sys.variance.resize(g_howManyFilts);

            for (size_t f = 0; f < g_howManyFilts; ++f) {
                sys.obsPhot[f] = obsDist(gen);
                sys.variance[f] = varDist(gen);
            }

            systems.push_back(sys);
        }

        // Prepare aligned arrays for SIMD (both x86 SSE and ARM NEON)
        size_t totalSize = nSystems * g_howManyFiltsAligned;

        sysVars.resize(totalSize);
        sysVar2.resize(totalSize);
        sysObs.resize(totalSize);
        sysFluxLo.resize(totalSize);
        sysFluxHi.resize(totalSize);

        const double TWO_M_PI = 2 * M_PI;
        const double nSD = 3.0;
        size_t idx = 0;
        for (const auto& sys : systems) {
            for (size_t k = 0; k < g_howManyFiltsAligned; ++k, ++idx) {
                if (k < g_howManyFilts) {
                    double var = sys.variance[k] * g_clust.varScale;
                    sysVars[idx] = 1.0 / var;
                    sysVar2[idx] = std::log(TWO_M_PI * var);
                    sysObs[idx] = sys.obsPhot[k];
                    double nSDTerm = nSD * sqrt(var);
                    sysFluxLo[idx] = exp10((sys.obsPhot[k] - nSDTerm) / -2.5);
                    sysFluxHi[idx] = exp10((sys.obsPhot[k] + nSDTerm) / -2.5);
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

    void TearDown(const benchmark::State&) override {
        systems.clear();
        sysVars.clear();
        sysVar2.clear();
        sysObs.clear();
        sysFluxLo.clear();
        sysFluxHi.clear();
    }

    vector<StellarSystem> systems;
    vector<double> sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi;
    size_t nSystems = 0;
    base::utility::ThreadPool pool{1};
};

BENCHMARK_DEFINE_F(NoBinariesFixture, MargEvolveNoBinaries_Reference)(benchmark::State& state) {
    if (!g_initialized) return;

    for (auto _ : state) {
        auto result = margEvolveNoBinaries_Reference(
            g_clust, *g_evoModels, *g_isochrone,
            systems, false, 0);
        benchmark::DoNotOptimize(result);
    }

    state.SetItemsProcessed(state.iterations() * nSystems);
    state.counters["stars"] = nSystems;
    state.counters["filters"] = g_howManyFilts;
}

BENCHMARK_REGISTER_F(NoBinariesFixture, MargEvolveNoBinaries_Reference)
    ->Unit(benchmark::kMillisecond)
    ->Arg(1)
    ->Arg(100)
    ->Arg(500)
    ->Arg(1000)
    ->Arg(2000);

BENCHMARK_DEFINE_F(NoBinariesFixture, MargEvolveNoBinaries)(benchmark::State& state) {
    if (!g_initialized) return;

    for (auto _ : state) {
        auto result = margEvolveNoBinaries(
            g_clust, *g_evoModels, *g_isochrone, pool,
            sysVars.data(), sysVar2.data(), sysObs.data(),
            nSystems, g_howManyFiltsAligned, g_howManyFilts,
            false, 0);
        benchmark::DoNotOptimize(result);
    }

    state.SetItemsProcessed(state.iterations() * nSystems);
    state.counters["stars"] = nSystems;
    state.counters["filters"] = g_howManyFilts;
}

BENCHMARK_REGISTER_F(NoBinariesFixture, MargEvolveNoBinaries)
    ->Unit(benchmark::kMillisecond)
    ->Arg(1)
    ->Arg(100)
    ->Arg(500)
    ->Arg(1000)
    ->Arg(2000);

// Binary benchmark - uses the same fixture but calls margEvolveWithBinary
BENCHMARK_DEFINE_F(NoBinariesFixture, MargEvolveWithBinary)(benchmark::State& state) {
    if (!g_initialized) return;

    for (auto _ : state) {
        auto result = margEvolveWithBinary(
            g_clust, systems, *g_evoModels, *g_isochrone, pool,
            false, 0);
        benchmark::DoNotOptimize(result);
    }

    state.SetItemsProcessed(state.iterations() * nSystems);
    state.counters["stars"] = nSystems;
    state.counters["filters"] = g_howManyFilts;
}

BENCHMARK_REGISTER_F(NoBinariesFixture, MargEvolveWithBinary)
    ->Unit(benchmark::kMillisecond)
    ->Arg(1)
    ->Arg(100)
    ->Arg(500)
    ->Arg(1000)
    ->Arg(2000);

// SIMD-optimized binary benchmark
BENCHMARK_DEFINE_F(NoBinariesFixture, MargEvolveWithBinary_SIMD)(benchmark::State& state) {
    if (!g_initialized) return;

    for (auto _ : state) {
        auto result = margEvolveWithBinary(
            g_clust, *g_evoModels, *g_isochrone, pool,
            systems,
            sysVars.data(), sysVar2.data(), sysObs.data(),
            sysFluxLo.data(), sysFluxHi.data(),
            g_howManyFiltsAligned, g_howManyFilts,
            false, 0);
        benchmark::DoNotOptimize(result);
    }

    state.SetItemsProcessed(state.iterations() * nSystems);
    state.counters["stars"] = nSystems;
    state.counters["filters"] = g_howManyFilts;
}

BENCHMARK_REGISTER_F(NoBinariesFixture, MargEvolveWithBinary_SIMD)
    ->Unit(benchmark::kMillisecond)
    ->Arg(1)
    ->Arg(100)
    ->Arg(500)
    ->Arg(1000)
    ->Arg(2000);

int main(int argc, char** argv) {
    // Separate benchmark args from BASE-9 args
    // BASE-9 args: -c <config>
    // benchmark args: everything else starting with --benchmark

    std::vector<char*> base9Args;
    std::vector<char*> benchArgs;

    base9Args.push_back(argv[0]);
    benchArgs.push_back(argv[0]);

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.rfind("--benchmark", 0) == 0) {
            benchArgs.push_back(argv[i]);
        } else {
            base9Args.push_back(argv[i]);
        }
    }

    try {
        initializeFromSettings(base9Args.size(), base9Args.data());
    } catch (const std::exception& e) {
        std::cerr << "Failed to initialize: " << e.what() << "\n";
        std::cerr << "Usage: " << argv[0] << " [-c config.yaml] [--benchmark_* options]\n";
        return 1;
    }

    int benchArgc = benchArgs.size();
    ::benchmark::Initialize(&benchArgc, benchArgs.data());
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();

    delete g_evoModels;
    delete g_isochrone;

    return 0;
}
