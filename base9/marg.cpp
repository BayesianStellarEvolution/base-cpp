#include <array>
#include <bitset>
#include <vector>

#include <cassert>
#include <cstdio>
#include <cmath>
#include <cstdlib>

#include "Cluster.hpp"
#include "Star.hpp"

#include "densities.hpp"
#include "GlobalSettings.hpp"
#include "Matrix.hpp"
#include "Model.hpp"
#include "Utility.hpp"
#include "WhiteDwarf.hpp"

using std::array;
using std::mutex;
using std::vector;
using base::utility::ThreadPool;

// ASSUMPTIONS
//   1. margEvolveWithBinary is not adversely affected if the StellarSystem is modified
//   2. margEvolveWithBinary passes the StellarSystem with primary.mass already set appropriately
//   3. margEvolveWithBinary expects secondary.mass to be overwritten, possibly immediately
static void calcPost (const double dMass, const Cluster &clust, vector<StellarSystem> &systems, const Model &evoModels, const Isochrone &isochrone, vector<double> &post, const vector<vector<double>> &secondaryMags, const bool modIsParallax)
{
    // This struct is only used inside this function, so we don't want it escaping
    struct diffStruct
    {
        diffStruct(unsigned short filter, double magLower, double magUpper)
            : filter(filter), magLower(magLower), magUpper(magUpper)
        {;}
        ~diffStruct() {;}

        const unsigned short filter;
        const double magLower;
        const double magUpper;
    };

    struct starDiff
    {
        starDiff(unsigned index, vector<diffStruct> diffs)
            : index(index), diffs(diffs)
        {;}
        ~starDiff() {;}

        const unsigned index;
        const vector<diffStruct> diffs;
    };

    // Useful constants
    const auto nSystems   = systems.size();
    const auto obsSize    = systems.front().obsPhot.size();

    const double logdMass = log (dMass); // This is a pre-optimzation so we only have to call log once (instead of nEntries times)

    // Rename system.primary.mass to primaryMass for this function
    const double primaryMass = systems.front().primary.mass;
    const double dMassRatio  = log (isochrone.eeps[0].mass / primaryMass);
    const double logPrior    = clust.logPriorMass (primaryMass);

    // Get the mags based on the primary's mass set by margEvolve
    // We call deriveCombinedMags here rather than primary.getMags so that we get abs/distMod correction
    // This keeps models which haven't been converted to absolute mags from breaking
    systems.front().secondary.mass = 0.0;
    const auto primaryMags  = systems.front().primary.getMags (clust, evoModels, isochrone);

    /**** now see if any binary companions are OK ****/
    const double nSD = 3.0;           // num of st dev from obs that will contribute to likelihood

    vector<starDiff> diffs;

    // Pre-calculate diffLow and diffUp so they can be cached for the loop below
    // In its own block to contain obsFilt and obsSize
    {
        vector<diffStruct> singleDiffs; // Keeps track of diffStructs for every filter
        singleDiffs.reserve(obsSize);

        const auto combinedMags = systems.front().deriveCombinedMags (clust, evoModels, isochrone, primaryMags, secondaryMags.back(), modIsParallax);

        // expTerms is pulled out of the loop, as it's the same across every star.
        vector<double> expTerms;

        for (size_t f = 0; f < obsSize; ++f)
        {
            expTerms.push_back(exp10 ((combinedMags[f] / -2.5)));
        }

        for (size_t s = 0; s < nSystems; ++s)
        {
            singleDiffs.clear();

            auto &system = systems[s];

            for (size_t f = 0; f < obsSize; ++f)
            {
                if (system.variance[f] > EPS)
                {
                    // Optimizations(?)
                    double obsPhot = system.obsPhot[f];               // Reduces array dereference penalty
                    double nSDTerm = nSD * sqrt (system.variance[f]); // Array dereference + multiplication and sqrt call

                    const double diffLow = exp10 (((obsPhot - nSDTerm) / -2.5)) - expTerms[f];
                    const double diffUp  = exp10 (((obsPhot + nSDTerm) / -2.5)) - expTerms[f];

                    if (diffLow <= 0.0 || diffLow == diffUp) // log(-x) == NaN
                    {
                        // Instead of returning post (a half-finished calculation), return 0 to signify
                        // that no mass combinations would lead to a positive posterior density
                        singleDiffs.clear();
                        break;
                    }
                    else
                    {
                        const double magLower = -2.5 * __builtin_log10 (diffLow);
                        const double magUpper = -2.5 * __builtin_log10 (diffUp);

                        assert(!std::isnan(magLower)); // Even though we checked it above, we still check it here
                        // magUpper is allowed to be NaN, as it is checked below.

                        singleDiffs.emplace_back(f, magLower, magUpper); // Make a diffStruct and stick it in diffs
                    }
                }
            }

            // Secondary mass to 0.0 to ensure we only look at the primary for the combinedMags.
            system.primary.mass   = primaryMass;
            system.secondary.mass = 0.0;

            if ( ! singleDiffs.empty() )
            {
                double tmpLogPost = system.logPost (clust, evoModels, isochrone, logPrior, primaryMags, secondaryMags.back(), modIsParallax)
                                  + logdMass + dMassRatio;

                post[s] += __builtin_exp (tmpLogPost);

                diffs.emplace_back(s, singleDiffs);
            }
        }
    }

    // Evaluate the posterior at all reasonable points in the isochrone

    auto isoSize = isochrone.eeps.size();
    vector<double> combinedMags;

    for (size_t i = 0; i < isoSize - 1; ++i)
    {
        // More optimization. Keeps us from dereferencing this array ~4 times.
        double isoMass = isochrone.eeps[i].mass;

        if (isoMass > primaryMass)
            break;

        if (!combinedMags.empty())
            combinedMags.clear();

        for (const auto &s : diffs )
        {
            StellarSystem &system = systems[s.index];

            // All filters should be acceptable if we are to evaluate at a grid point
            // okMass starts out true, and gets set false if any filter is not acceptable
            bool okMass = true;

            // For all the diffSets (i.e., every filter)
            for (const auto &diff : s.diffs)
            {
                // The isochrone mag should be between magLower and magUpper
                // and the isochrone mass should be less than primaryMass
                //
                // Exception: If magUpper is NaN (i.e., it was < 0 as diffUp above), consider it irrelevant
                double mag = isochrone.eeps[i].mags[diff.filter];

                if ( ! (mag >= diff.magLower
                        && (mag <= diff.magUpper || std::isnan(diff.magUpper))))
                {
                    okMass = false; // If it isn't, that isn't OK
                    break;          // And break out of the loop to avoid further calculation (pre-optimization)
                }
            }

            // Now, only if we left the above loop without triggering the condition, calculate the posterior density for the system
            if (okMass)
            {
                system.setMassRatio (isoMass / primaryMass); // Set the massRatio (and therefore the secondary mass)

                if (combinedMags.empty())
                    combinedMags = system.deriveCombinedMags (clust, evoModels, isochrone, primaryMags, secondaryMags[i], modIsParallax);

                // Calculate the posterior density for this system
                double tmpLogPost = system.logPost (clust, evoModels, isochrone, logPrior, combinedMags)
                    + logdMass
                    + __builtin_log ((isochrone.eeps[i + 1].mass - isoMass) / primaryMass);

                post[s.index] += __builtin_exp (tmpLogPost); // And add it to the running total
            }
        }
    }
}


/* evaluate on a grid of primary mass and mass ratio to approximate the integral */
vector<double> margEvolveWithBinary (const Cluster &clust, vector<StellarSystem> &systems, const Model &evoModels, const Isochrone &isochrone, ThreadPool &, const bool modIsParallax, const int eepInterpolationPower)
{
    assert(isochrone.eeps.size() >= 2);

    vector<vector<double>> secondaryMags;

    for (const auto &e : isochrone.eeps)
    {
        Star s;
        s.mass = e.mass;

        secondaryMags.push_back(s.getMags(clust, evoModels, isochrone));
    }

    // secondaryMags.back() is the 0-mass magnitudes
    {
        Star s;
        s.mass = 0.0;

        secondaryMags.push_back(s.getMags(clust, evoModels, isochrone));
    }

    const int isoIncrem = GlobalSettings::calcEepSteps(eepInterpolationPower);

    vector<double> post(systems.size(), 0.0);

    {
        auto isoSize = isochrone.eeps.size();

        for (size_t m = 0; m < isoSize - 1; ++m)
        {
            double dIsoMass = isochrone.eeps[m + 1].mass - isochrone.eeps[m].mass;

            // In the event that we have an invalid range, skip that range
            // This generally occurs only at very high EEPs, where the masses are close together
            if (dIsoMass < 0.0) continue;

            double dMass    = dIsoMass / isoIncrem;

            for (auto k = 0; k < isoIncrem; ++k)
            {
                systems.front().primary.mass = isochrone.eeps[m].mass + (k * dMass);

                calcPost (dMass, clust, systems, evoModels, isochrone, post, secondaryMags, modIsParallax);
            }
        }
    }

    return post;
}

#ifdef IS_ARM
#include <arm_neon.h>
#define ROUND_DOWN(x, s) ((x) & ~((s)-1))

// SIMD-optimized likelihood calculation for a single system
// Uses pre-computed 1/var, log(2πvar), and observation arrays
inline static double calcLogLikelihood_neon(
    const double* sVar,      // 1.0 / variance
    const double* sVa2,      // log(2π * variance)
    const double* sObs,      // observed photometry
    const double* mags,      // model magnitudes (aligned)
    const size_t obsUnaligned)
{
    const float64x2_t half = vdupq_n_f64(-0.5);
    float64x2_t filterSum = vdupq_n_f64(0.0);

    size_t f;
    for (f = 0; f < ROUND_DOWN(obsUnaligned, 2); f += 2)
    {
        const float64x2_t magF = vld1q_f64(mags + f);
        const float64x2_t varF = vld1q_f64(sVar + f);
        const float64x2_t va2F = vld1q_f64(sVa2 + f);
        const float64x2_t obsF = vld1q_f64(sObs + f);

        float64x2_t deltaMags = vsubq_f64(magF, obsF);
        deltaMags = vmulq_f64(deltaMags, deltaMags);
        deltaMags = vmulq_f64(deltaMags, varF);
        deltaMags = vaddq_f64(deltaMags, va2F);
        deltaMags = vmulq_f64(deltaMags, half);

        filterSum = vaddq_f64(filterSum, deltaMags);
    }

    double lik = vgetq_lane_f64(filterSum, 0) + vgetq_lane_f64(filterSum, 1);

    // Handle remaining odd filter
    for (; f < obsUnaligned; ++f)
    {
        double deltaMags = mags[f] - sObs[f];
        deltaMags = deltaMags * deltaMags * sVar[f] + sVa2[f];
        lik += deltaMags * (-0.5);
    }

    return lik;
}

// SIMD-optimized calcPost with pre-computed terms
static void calcPost_simd(
    const double dMass,
    const Cluster &clust,
    vector<StellarSystem> &systems,
    const Model &evoModels,
    const Isochrone &isochrone,
    vector<double> &post,
    const vector<vector<double>> &secondaryMags,
    const vector<vector<double>> &secondaryFluxes,  // Pre-computed exp10(secondaryMags[i][f] / -2.5)
    const double* const distCorrection,              // Pre-computed distance + (absCoeffs[f] - 1) * abs
    const double* const systemVars,
    const double* const systemVar2,
    const double* const systemObs,
    const double* const systemFluxLo,   // Pre-computed exp10((obs - 3*sqrt(var)) / -2.5)
    const double* const systemFluxHi,   // Pre-computed exp10((obs + 3*sqrt(var)) / -2.5)
    const size_t obsSize,
    const size_t obsUnaligned)
{
    struct diffStruct
    {
        diffStruct(unsigned short filter, double magLower, double magUpper)
            : filter(filter), magLower(magLower), magUpper(magUpper) {}
        const unsigned short filter;
        const double magLower;
        const double magUpper;
    };

    struct starDiff
    {
        starDiff(unsigned index, vector<diffStruct> diffs)
            : index(index), diffs(diffs) {}
        const unsigned index;
        const vector<diffStruct> diffs;
    };

    const auto nSystems = systems.size();
    const double logdMass = __builtin_log(dMass);
    const double primaryMass = systems.front().primary.mass;
    const double dMassRatio = __builtin_log(isochrone.eeps[0].mass / primaryMass);
    const double logPrior = clust.logPriorMass(primaryMass);

    // Get primary magnitudes (from isochrone interpolation)
    Star pStar;
    pStar.mass = primaryMass;
    const auto primaryMags = pStar.msRgbEvol(isochrone);

    // Pre-compute primary fluxes
    vector<double> primaryFluxes(obsUnaligned);
    for (size_t f = 0; f < obsUnaligned; ++f)
        primaryFluxes[f] = exp10(primaryMags[f] / -2.5);

    const double nSD = 3.0;
    vector<starDiff> diffs;

    // Aligned buffer for combined magnitudes
    double* aligned_mags = reinterpret_cast<double*>(aligned_alloc(16, obsSize * sizeof(double)));

    // q=0 case: combined mags = primary mags + distance correction (no secondary contribution)
    {
        // Compute combined mags for q=0 (secondary flux = 0, use secondaryFluxes.back() which is for mass=0)
        for (size_t f = 0; f < obsUnaligned; ++f)
        {
            double flux = primaryFluxes[f] + secondaryFluxes.back()[f];
            aligned_mags[f] = -2.5 * __builtin_log10(flux) + distCorrection[f];
        }

        // Compute expTerms for filtering
        vector<double> expTerms(obsUnaligned);
        for (size_t f = 0; f < obsUnaligned; ++f)
            expTerms[f] = exp10(aligned_mags[f] / -2.5);

        vector<diffStruct> singleDiffs;
        singleDiffs.reserve(obsUnaligned);

        for (size_t s = 0; s < nSystems; ++s)
        {
            singleDiffs.clear();
            auto &system = systems[s];
            const double* sFluxLo = systemFluxLo + s * obsSize;
            const double* sFluxHi = systemFluxHi + s * obsSize;

            for (size_t f = 0; f < obsUnaligned; ++f)
            {
                if (system.variance[f] > EPS)
                {
                    const double diffLow = sFluxLo[f] - expTerms[f];
                    const double diffUp = sFluxHi[f] - expTerms[f];

                    if (diffLow <= 0.0 || diffLow == diffUp)
                    {
                        singleDiffs.clear();
                        break;
                    }
                    else
                    {
                        const double magLower = -2.5 * __builtin_log10(diffLow);
                        const double magUpper = -2.5 * __builtin_log10(diffUp);
                        singleDiffs.emplace_back(f, magLower, magUpper);
                    }
                }
            }

            system.primary.mass = primaryMass;
            system.secondary.mass = 0.0;

            if (!singleDiffs.empty())
            {
                const double* sVar = systemVars + s * obsSize;
                const double* sVa2 = systemVar2 + s * obsSize;
                const double* sObs = systemObs + s * obsSize;

                double likelihood = calcLogLikelihood_neon(sVar, sVa2, sObs, aligned_mags, obsUnaligned);
                double tmpLogPost = logPrior + likelihood + logdMass + dMassRatio;

                post[s] += __builtin_exp(tmpLogPost);
                diffs.emplace_back(s, singleDiffs);
            }
        }
    }

    // Non-zero mass ratio cases
    auto isoSize = isochrone.eeps.size();

    for (size_t i = 0; i < isoSize - 1; ++i)
    {
        double isoMass = isochrone.eeps[i].mass;
        if (isoMass > primaryMass)
            break;

        bool combinedMagsComputed = false;

        for (const auto &s : diffs)
        {
            bool okMass = true;

            for (const auto &diff : s.diffs)
            {
                double mag = isochrone.eeps[i].mags[diff.filter];
                if (!(mag >= diff.magLower && (mag <= diff.magUpper || std::isnan(diff.magUpper))))
                {
                    okMass = false;
                    break;
                }
            }

            if (okMass)
            {
                // Compute combined mags once per isochrone point (shared across valid systems)
                if (!combinedMagsComputed)
                {
                    for (size_t f = 0; f < obsUnaligned; ++f)
                    {
                        double flux = primaryFluxes[f] + secondaryFluxes[i][f];
                        aligned_mags[f] = -2.5 * __builtin_log10(flux) + distCorrection[f];
                    }
                    combinedMagsComputed = true;
                }

                const double* sVar = systemVars + s.index * obsSize;
                const double* sVa2 = systemVar2 + s.index * obsSize;
                const double* sObs = systemObs + s.index * obsSize;

                double likelihood = calcLogLikelihood_neon(sVar, sVa2, sObs, aligned_mags, obsUnaligned);
                double tmpLogPost = logPrior + likelihood + logdMass
                    + __builtin_log((isochrone.eeps[i + 1].mass - isoMass) / primaryMass);

                post[s.index] += __builtin_exp(tmpLogPost);
            }
        }
    }

    free(aligned_mags);
}

// SIMD-optimized margEvolveWithBinary
vector<double> margEvolveWithBinary(
    const Cluster &clust,
    const Model &evoModels,
    const Isochrone &isochrone,
    ThreadPool &,
    vector<StellarSystem> &systems,
    const double* const systemVars,
    const double* const systemVar2,
    const double* const systemObs,
    const double* const systemFluxLo,
    const double* const systemFluxHi,
    const size_t obsSize,
    const size_t obsUnaligned,
    const bool modIsParallax,
    const int eepInterpolationPower)
{
    assert(isochrone.eeps.size() >= 2);

    // Pre-compute secondary magnitudes for all isochrone points
    vector<vector<double>> secondaryMags;
    for (const auto &e : isochrone.eeps)
    {
        Star s;
        s.mass = e.mass;
        secondaryMags.push_back(s.getMags(clust, evoModels, isochrone));
    }
    {
        Star s;
        s.mass = 0.0;
        secondaryMags.push_back(s.getMags(clust, evoModels, isochrone));
    }

    // Pre-compute secondary fluxes: exp10(secondaryMags[i][f] / -2.5)
    vector<vector<double>> secondaryFluxes;
    for (const auto& mags : secondaryMags)
    {
        vector<double> fluxes(obsUnaligned);
        for (size_t f = 0; f < obsUnaligned; ++f)
            fluxes[f] = exp10(mags[f] / -2.5);
        secondaryFluxes.push_back(std::move(fluxes));
    }

    // Pre-compute distance correction: distance + (absCoeffs[f] - 1) * abs
    double distance;
    if (modIsParallax)
    {
        auto parsecs = 1.0 / clust.mod;
        distance = clust.abs + 5.0 * log10(parsecs) - 5.0;
    }
    else
    {
        distance = clust.mod;
    }

    double* distCorrection = reinterpret_cast<double*>(aligned_alloc(16, obsSize * sizeof(double)));
    for (size_t f = 0; f < obsUnaligned; ++f)
        distCorrection[f] = distance + (evoModels.absCoeffs[f] - 1.0) * clust.abs;

    const int isoIncrem = GlobalSettings::calcEepSteps(eepInterpolationPower);

    vector<double> post(systems.size(), 0.0);

    {
        auto isoSize = isochrone.eeps.size();

        for (size_t m = 0; m < isoSize - 1; ++m)
        {
            double dIsoMass = isochrone.eeps[m + 1].mass - isochrone.eeps[m].mass;

            if (dIsoMass < 0.0) continue;

            double dMass = dIsoMass / isoIncrem;

            for (auto k = 0; k < isoIncrem; ++k)
            {
                systems.front().primary.mass = isochrone.eeps[m].mass + (k * dMass);

                calcPost_simd(dMass, clust, systems, evoModels, isochrone, post, secondaryMags,
                             secondaryFluxes, distCorrection,
                             systemVars, systemVar2, systemObs, systemFluxLo, systemFluxHi,
                             obsSize, obsUnaligned);
            }
        }
    }

    free(distCorrection);
    return post;
}

#else
// x86 SSE version
#include <emmintrin.h>
#include <pmmintrin.h>
#define ROUND_DOWN_MARG(x, s) ((x) & ~((s)-1))

inline static double calcLogLikelihood_sse(
    const double* sVar,
    const double* sVa2,
    const double* sObs,
    const double* mags,
    const size_t obsUnaligned)
{
    const __m128d half = _mm_set1_pd(-0.5);
    __m128d filterSum = _mm_setzero_pd();

    size_t f;
    for (f = 0; f < ROUND_DOWN_MARG(obsUnaligned, 2); f += 2)
    {
        const __m128d magF = _mm_load_pd(mags + f);
        const __m128d varF = _mm_load_pd(sVar + f);
        const __m128d va2F = _mm_load_pd(sVa2 + f);
        const __m128d obsF = _mm_load_pd(sObs + f);

        __m128d deltaMags = _mm_sub_pd(magF, obsF);
        deltaMags = _mm_mul_pd(deltaMags, deltaMags);
        deltaMags = _mm_mul_pd(deltaMags, varF);
        deltaMags = _mm_add_pd(deltaMags, va2F);
        deltaMags = _mm_mul_pd(deltaMags, half);

        filterSum = _mm_add_pd(filterSum, deltaMags);
    }

    filterSum = _mm_hadd_pd(filterSum, filterSum);
    double lik;
    _mm_store_sd(&lik, filterSum);

    for (; f < obsUnaligned; ++f)
    {
        double deltaMags = mags[f] - sObs[f];
        deltaMags = deltaMags * deltaMags * sVar[f] + sVa2[f];
        lik += deltaMags * (-0.5);
    }

    return lik;
}

static void calcPost_simd(
    const double dMass,
    const Cluster &clust,
    vector<StellarSystem> &systems,
    const Model &evoModels,
    const Isochrone &isochrone,
    vector<double> &post,
    const vector<vector<double>> &secondaryMags,
    const vector<vector<double>> &secondaryFluxes,
    const double* const distCorrection,
    const double* const systemVars,
    const double* const systemVar2,
    const double* const systemObs,
    const double* const systemFluxLo,
    const double* const systemFluxHi,
    const size_t obsSize,
    const size_t obsUnaligned)
{
    struct diffStruct
    {
        diffStruct(unsigned short filter, double magLower, double magUpper)
            : filter(filter), magLower(magLower), magUpper(magUpper) {}
        const unsigned short filter;
        const double magLower;
        const double magUpper;
    };

    struct starDiff
    {
        starDiff(unsigned index, vector<diffStruct> diffs)
            : index(index), diffs(diffs) {}
        const unsigned index;
        const vector<diffStruct> diffs;
    };

    const auto nSystems = systems.size();
    const double logdMass = __builtin_log(dMass);
    const double primaryMass = systems.front().primary.mass;
    const double dMassRatio = __builtin_log(isochrone.eeps[0].mass / primaryMass);
    const double logPrior = clust.logPriorMass(primaryMass);

    // Get primary magnitudes (from isochrone interpolation)
    Star pStar;
    pStar.mass = primaryMass;
    const auto primaryMags = pStar.msRgbEvol(isochrone);

    // Pre-compute primary fluxes
    vector<double> primaryFluxes(obsUnaligned);
    for (size_t f = 0; f < obsUnaligned; ++f)
        primaryFluxes[f] = exp10(primaryMags[f] / -2.5);

    const double nSD = 3.0;
    vector<starDiff> diffs;

    double* aligned_mags = reinterpret_cast<double*>(_mm_malloc(obsSize * sizeof(double), 16));

    // q=0 case
    {
        for (size_t f = 0; f < obsUnaligned; ++f)
        {
            double flux = primaryFluxes[f] + secondaryFluxes.back()[f];
            aligned_mags[f] = -2.5 * __builtin_log10(flux) + distCorrection[f];
        }

        vector<double> expTerms(obsUnaligned);
        for (size_t f = 0; f < obsUnaligned; ++f)
            expTerms[f] = exp10(aligned_mags[f] / -2.5);

        vector<diffStruct> singleDiffs;
        singleDiffs.reserve(obsUnaligned);

        for (size_t s = 0; s < nSystems; ++s)
        {
            singleDiffs.clear();
            auto &system = systems[s];
            const double* sFluxLo = systemFluxLo + s * obsSize;
            const double* sFluxHi = systemFluxHi + s * obsSize;

            for (size_t f = 0; f < obsUnaligned; ++f)
            {
                if (system.variance[f] > EPS)
                {
                    const double diffLow = sFluxLo[f] - expTerms[f];
                    const double diffUp = sFluxHi[f] - expTerms[f];

                    if (diffLow <= 0.0 || diffLow == diffUp)
                    {
                        singleDiffs.clear();
                        break;
                    }
                    else
                    {
                        const double magLower = -2.5 * __builtin_log10(diffLow);
                        const double magUpper = -2.5 * __builtin_log10(diffUp);
                        singleDiffs.emplace_back(f, magLower, magUpper);
                    }
                }
            }

            system.primary.mass = primaryMass;
            system.secondary.mass = 0.0;

            if (!singleDiffs.empty())
            {
                const double* sVar = systemVars + s * obsSize;
                const double* sVa2 = systemVar2 + s * obsSize;
                const double* sObs = systemObs + s * obsSize;

                double likelihood = calcLogLikelihood_sse(sVar, sVa2, sObs, aligned_mags, obsUnaligned);
                double tmpLogPost = logPrior + likelihood + logdMass + dMassRatio;

                post[s] += __builtin_exp(tmpLogPost);
                diffs.emplace_back(s, singleDiffs);
            }
        }
    }

    auto isoSize = isochrone.eeps.size();

    for (size_t i = 0; i < isoSize - 1; ++i)
    {
        double isoMass = isochrone.eeps[i].mass;
        if (isoMass > primaryMass)
            break;

        bool combinedMagsComputed = false;

        for (const auto &s : diffs)
        {
            bool okMass = true;

            for (const auto &diff : s.diffs)
            {
                double mag = isochrone.eeps[i].mags[diff.filter];
                if (!(mag >= diff.magLower && (mag <= diff.magUpper || std::isnan(diff.magUpper))))
                {
                    okMass = false;
                    break;
                }
            }

            if (okMass)
            {
                if (!combinedMagsComputed)
                {
                    for (size_t f = 0; f < obsUnaligned; ++f)
                    {
                        double flux = primaryFluxes[f] + secondaryFluxes[i][f];
                        aligned_mags[f] = -2.5 * __builtin_log10(flux) + distCorrection[f];
                    }
                    combinedMagsComputed = true;
                }

                const double* sVar = systemVars + s.index * obsSize;
                const double* sVa2 = systemVar2 + s.index * obsSize;
                const double* sObs = systemObs + s.index * obsSize;

                double likelihood = calcLogLikelihood_sse(sVar, sVa2, sObs, aligned_mags, obsUnaligned);
                double tmpLogPost = logPrior + likelihood + logdMass
                    + __builtin_log((isochrone.eeps[i + 1].mass - isoMass) / primaryMass);

                post[s.index] += __builtin_exp(tmpLogPost);
            }
        }
    }

    _mm_free(aligned_mags);
}

vector<double> margEvolveWithBinary(
    const Cluster &clust,
    const Model &evoModels,
    const Isochrone &isochrone,
    ThreadPool &,
    vector<StellarSystem> &systems,
    const double* const systemVars,
    const double* const systemVar2,
    const double* const systemObs,
    const double* const systemFluxLo,
    const double* const systemFluxHi,
    const size_t obsSize,
    const size_t obsUnaligned,
    const bool modIsParallax,
    const int eepInterpolationPower)
{
    assert(isochrone.eeps.size() >= 2);

    // Pre-compute secondary magnitudes for all isochrone points
    vector<vector<double>> secondaryMags;
    for (const auto &e : isochrone.eeps)
    {
        Star s;
        s.mass = e.mass;
        secondaryMags.push_back(s.getMags(clust, evoModels, isochrone));
    }
    {
        Star s;
        s.mass = 0.0;
        secondaryMags.push_back(s.getMags(clust, evoModels, isochrone));
    }

    // Pre-compute secondary fluxes: exp10(secondaryMags[i][f] / -2.5)
    vector<vector<double>> secondaryFluxes;
    for (const auto& mags : secondaryMags)
    {
        vector<double> fluxes(obsUnaligned);
        for (size_t f = 0; f < obsUnaligned; ++f)
            fluxes[f] = exp10(mags[f] / -2.5);
        secondaryFluxes.push_back(std::move(fluxes));
    }

    // Pre-compute distance correction: distance + (absCoeffs[f] - 1) * abs
    double distance;
    if (modIsParallax)
    {
        auto parsecs = 1.0 / clust.mod;
        distance = clust.abs + 5.0 * log10(parsecs) - 5.0;
    }
    else
    {
        distance = clust.mod;
    }

    double* distCorrection = reinterpret_cast<double*>(_mm_malloc(obsSize * sizeof(double), 16));
    for (size_t f = 0; f < obsUnaligned; ++f)
        distCorrection[f] = distance + (evoModels.absCoeffs[f] - 1.0) * clust.abs;

    const int isoIncrem = GlobalSettings::calcEepSteps(eepInterpolationPower);

    vector<double> post(systems.size(), 0.0);

    {
        auto isoSize = isochrone.eeps.size();

        for (size_t m = 0; m < isoSize - 1; ++m)
        {
            double dIsoMass = isochrone.eeps[m + 1].mass - isochrone.eeps[m].mass;

            if (dIsoMass < 0.0) continue;

            double dMass = dIsoMass / isoIncrem;

            for (auto k = 0; k < isoIncrem; ++k)
            {
                systems.front().primary.mass = isochrone.eeps[m].mass + (k * dMass);

                calcPost_simd(dMass, clust, systems, evoModels, isochrone, post, secondaryMags,
                             secondaryFluxes, distCorrection,
                             systemVars, systemVar2, systemObs, systemFluxLo, systemFluxHi,
                             obsSize, obsUnaligned);
            }
        }
    }

    _mm_free(distCorrection);
    return post;
}
#endif
