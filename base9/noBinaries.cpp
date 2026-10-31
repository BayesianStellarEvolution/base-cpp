#include <algorithm>
#include <array>
#include <bitset>
#include <iostream>
#include <vector>

#include <cassert>
#include <cstdio>
#include <cmath>
#include <cstdlib>

#include "Cluster.hpp"
#include "LinearTransform.hpp"
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

// Forward declaration of scalar reference (defined at end of file)
vector<double> margEvolveNoBinaries_Reference(const Cluster &, const Model &, const Isochrone &, vector<StellarSystem> &, const bool, const int);

#ifndef IS_ARM

#include <emmintrin.h> // SSE2 (_pd, _m128d)
#include <pmmintrin.h> // SSE3 (hadd)
#include <sleef.h>

#define ROUND_DOWN(x, s) ((x) & ~((s)-1))

inline static __m128d calcDeltaMags(const double* sVar, const double* sVa2, const double* sObs, const __m128d magF, const __m128d half, const size_t f)
{
    const __m128d varF = _mm_load_pd(sVar + f);
    const __m128d va2F = _mm_load_pd(sVa2 + f);
    const __m128d obsF = _mm_load_pd(sObs + f);

    // Do the math
    __m128d deltaMags = _mm_sub_pd(magF, obsF);   // magsF - obsF
    deltaMags = _mm_mul_pd(deltaMags, deltaMags); // sqr (deltaMags)
    deltaMags = _mm_mul_pd(deltaMags, varF);      // deltaMags * varF
    deltaMags = _mm_add_pd(deltaMags, va2F);      // deltaMags + va2F
    deltaMags = _mm_mul_pd(deltaMags, half);      // deltaMags * (-0.5)

    return deltaMags;
}

inline static __m128d calcSingleFilterResult(const double* sVar1, const double* sVa21, const double* sObs1, const double* sVar2, const double* sVa22, const double* sObs2, const __m128d magF, const __m128d half, const size_t f)
{
    // Load the low element with star 1 in the low bits
    __m128d varF = _mm_load_sd(sVar1 + f);
    __m128d va2F = _mm_load_sd(sVa21 + f);
    __m128d obsF = _mm_load_sd(sObs1 + f);

    // And with star 2 in the high bits
    varF = _mm_loadh_pd(varF, sVar2 + f);
    va2F = _mm_loadh_pd(va2F, sVa22 + f);
    obsF = _mm_loadh_pd(obsF, sObs2 + f);

    // Do the normal math for one star
    __m128d deltaMags = _mm_sub_pd(magF, obsF);   // magsF - obsF
    deltaMags = _mm_mul_pd(deltaMags, deltaMags); // sqr (deltaMags)
    deltaMags = _mm_mul_pd(deltaMags, varF);      // deltaMags * varF
    deltaMags = _mm_add_pd(deltaMags, va2F);      // deltaMags + va2F
    deltaMags = _mm_mul_pd(deltaMags, half);      // deltaMags * (-0.5)

    return deltaMags;
}

inline static void addLogLikelihood(const __m128d likelihood, double* sPost)
{
    const __m128d posts = Sleef_expd2_u10sse2(likelihood);

    __m128d tPost = _mm_load_pd(sPost);

    tPost = _mm_add_pd(tPost, posts);
    _mm_store_pd(sPost, tPost);
}


// Memory layout
//   struct threadStorage {
//     double combinedMags[obsSize]
//     double systems[nSystems]
//   };
//
//   threadStorage storage[nThreads]
//   double post[nSystems]
vector<double> margEvolveNoBinaries(const Cluster &clust, const Model &evoModels, const Isochrone &isochrone, ThreadPool &, const double* const systemVars, const double* const systemVar2, const double* const systemObs, const size_t nSystems, const size_t obsSize, const size_t obsUnaligned, const bool modIsParallax, const int eepInterpolationPower)
{
    const int isoIncrem = GlobalSettings::calcEepSteps(eepInterpolationPower);

    double* post = reinterpret_cast<double*>(_mm_malloc(nSystems * sizeof(double), 16));

    // Correct distance if it's in parallax
    double distance;

    if (modIsParallax)
    {
        auto parsecs = 1 / clust.mod;
        distance = clust.abs + 5 * log10(parsecs) - 5;
    }
    else
    {
        distance = clust.mod;
    }

    {
        size_t i;
        auto tPost = post;
        const __m128d zero = _mm_setzero_pd();

        for (i = 0; i < ROUND_DOWN(nSystems, 8); i += 8, tPost += 8)
        {
            _mm_stream_pd(tPost,     zero);
            _mm_stream_pd(tPost + 2, zero);
            _mm_stream_pd(tPost + 4, zero);
            _mm_stream_pd(tPost + 6, zero);
        }
        for (; i < nSystems; ++i)
            post[i] = 0;
    }

    double* aligned_mags = reinterpret_cast<double*>(_mm_malloc(obsSize * sizeof(double), 16));

    {
        auto agbTipMass = isochrone.agbTipMass();
        auto isoSize    = isochrone.eeps.size();

        for (size_t m = 0; m < isoSize - 1; ++m)
        {
            double dIsoMass = isochrone.eeps[m + 1].mass - isochrone.eeps[m].mass;

            // In the event that we have an invalid range, skip that range
            // This generally occurs only at very high EEPs, where the masses are close together
            if (dIsoMass > 0.0)
            {
                Star s;
                vector<double> combinedMags;

                double dMass    = dIsoMass / isoIncrem;
                double logdMass = __builtin_log (dMass);

                double cMass = isochrone.eeps[m].mass;

                for (auto k = 0; k < isoIncrem; ++k)
                {
                    double primaryMass = cMass + (k * dMass);

                    if (primaryMass >= agbTipMass)
                    {
                        primaryMass = agbTipMass;

                        // Kludgey way of exiting the for loop
                        // Only works if we don't use K below here
                        // Ensures we don't add the agbTip likelihood more than once
                        k = isoIncrem;
                    }

                    const double logPrior = clust.logPriorMass (primaryMass);

                    s.mass       = primaryMass;
                    combinedMags = s.msRgbEvol(isochrone);

                    for (size_t f = 0; f < obsUnaligned; ++f)
                    {
                        combinedMags[f] += distance;
                        combinedMags[f] += (evoModels.absCoeffs[f] - 1.0) * clust.abs;       // add A_[u-k] (standard defn of modulus already includes Av)
                    }

                    std::copy(combinedMags.begin(), combinedMags.end(), aligned_mags);

                    // The first star in the loop
                    auto sVar1 = systemVars;
                    auto sVa21 = systemVar2;
                    auto sObs1 = systemObs;

                    auto sVar2 = sVar1 + obsSize;
                    auto sVa22 = sVa21 + obsSize;
                    auto sObs2 = sObs1 + obsSize;

                    auto sVar3 = sVar2 + obsSize;
                    auto sVa23 = sVa22 + obsSize;
                    auto sObs3 = sObs2 + obsSize;

                    auto sVar4 = sVar3 + obsSize;
                    auto sVa24 = sVa23 + obsSize;
                    auto sObs4 = sObs3 + obsSize;

                    auto sVar5 = sVar4 + obsSize;
                    auto sVa25 = sVa24 + obsSize;
                    auto sObs5 = sObs4 + obsSize;

                    auto sVar6 = sVar5 + obsSize;
                    auto sVa26 = sVa25 + obsSize;
                    auto sObs6 = sObs5 + obsSize;

                    auto sVar7 = sVar6 + obsSize;
                    auto sVa27 = sVa26 + obsSize;
                    auto sObs7 = sObs6 + obsSize;

                    auto sVar8 = sVar7 + obsSize;
                    auto sVa28 = sVa27 + obsSize;
                    auto sObs8 = sObs7 + obsSize;

                    auto sPost = post;

                    size_t s;

                    // We calculate 4 stars at a time, so it has to 4x twice the number of filters
                    const auto perIter = 8;
                    const auto skip = obsSize * perIter; // % 2 to add an extra and make up for padding

                    const __m128d half   = _mm_set1_pd(-0.5);
                    const __m128d priors = _mm_set1_pd(logPrior + logdMass);

                    for (s = 0; s < ROUND_DOWN(nSystems, perIter); s += perIter, sPost += perIter,
                             sVar1 += skip, sVa21 += skip, sObs1 += skip,
                             sVar2 += skip, sVa22 += skip, sObs2 += skip,
                             sVar3 += skip, sVa23 += skip, sObs3 += skip,
                             sVar4 += skip, sVa24 += skip, sObs4 += skip,
                             sVar5 += skip, sVa25 += skip, sObs5 += skip,
                             sVar6 += skip, sVa26 += skip, sObs6 += skip,
                             sVar7 += skip, sVa27 += skip, sObs7 += skip,
                             sVar8 += skip, sVa28 += skip, sObs8 += skip)
                    {
                        __m128d likelihood1 = priors;
                        __m128d likelihood2 = priors;
                        __m128d likelihood3 = priors;
                        __m128d likelihood4 = priors;

                        size_t f;

                        // Main fun for filters. 4 stars, 2 filters per loop
                        for (f = 0; f < ROUND_DOWN(obsUnaligned, 2); f += 2)
                        {
                            // This is shared
                            const __m128d magF = _mm_load_pd(aligned_mags + f);

                            // Do the math
                            __m128d deltaMags1 = calcDeltaMags(sVar1, sVa21, sObs1, magF, half, f);
                            __m128d deltaMags2 = calcDeltaMags(sVar2, sVa22, sObs2, magF, half, f);
                            __m128d deltaMags3 = calcDeltaMags(sVar3, sVa23, sObs3, magF, half, f);
                            __m128d deltaMags4 = calcDeltaMags(sVar4, sVa24, sObs4, magF, half, f);
                            __m128d deltaMags5 = calcDeltaMags(sVar5, sVa25, sObs5, magF, half, f);
                            __m128d deltaMags6 = calcDeltaMags(sVar6, sVa26, sObs6, magF, half, f);
                            __m128d deltaMags7 = calcDeltaMags(sVar7, sVa27, sObs7, magF, half, f);
                            __m128d deltaMags8 = calcDeltaMags(sVar8, sVa28, sObs8, magF, half, f);

                            // Horizontally add the two stars worth of filters
                            const __m128d result1 = _mm_hadd_pd(deltaMags1, deltaMags2);
                            const __m128d result2 = _mm_hadd_pd(deltaMags3, deltaMags4);
                            const __m128d result3 = _mm_hadd_pd(deltaMags5, deltaMags6);
                            const __m128d result4 = _mm_hadd_pd(deltaMags7, deltaMags8);

                            // Sum the likelihoods
                            likelihood1 = _mm_add_pd(likelihood1, result1);
                            likelihood2 = _mm_add_pd(likelihood2, result2);
                            likelihood3 = _mm_add_pd(likelihood3, result3);
                            likelihood4 = _mm_add_pd(likelihood4, result4);
                        }

                        // Clean up the last filter in 2 stars (for odd numbers of filters)
                        for (; f < obsUnaligned; ++f)
                        {
                            // Broadcast magF
                            const __m128d magF = _mm_load1_pd(aligned_mags + f);

                            const __m128d result1 = calcSingleFilterResult(sVar1, sVa21, sObs1, sVar2, sVa22, sObs2, magF, half, f);
                            const __m128d result2 = calcSingleFilterResult(sVar3, sVa23, sObs3, sVar4, sVa24, sObs4, magF, half, f);
                            const __m128d result3 = calcSingleFilterResult(sVar5, sVa25, sObs5, sVar6, sVa26, sObs6, magF, half, f);
                            const __m128d result4 = calcSingleFilterResult(sVar7, sVa27, sObs7, sVar8, sVa28, sObs8, magF, half, f);

                            // Pretend we did the hadd and sum with likelihood
                            likelihood1 = _mm_add_pd(likelihood1, result1);
                            likelihood2 = _mm_add_pd(likelihood2, result2);
                            likelihood3 = _mm_add_pd(likelihood3, result3);
                            likelihood4 = _mm_add_pd(likelihood4, result4);
                        }


                        addLogLikelihood(likelihood1, sPost);
                        addLogLikelihood(likelihood2, sPost + 2);
                        addLogLikelihood(likelihood3, sPost + 4);
                        addLogLikelihood(likelihood4, sPost + 6);
                    }

                    // Clean up the last star
                    // This is automatically turned into SSE on non-AVX platforms
                    // Arguably, it shouldn't be SSE on AVX platforms unless we _mm_zeroupper
                    for (; s < nSystems; ++s, ++sPost,
                             sVar1 += obsSize, sVa21 += obsSize, sObs1 += obsSize)
                    {
                        __m128d likelihood = priors;

                        for (size_t f = 0; f < obsUnaligned; ++f)
                        {
                            const __m128d varF = _mm_load_sd(sVar1 + f);
                            const __m128d va2F = _mm_load_sd(sVa21 + f);
                            const __m128d obsF = _mm_load_sd(sObs1 + f);
                            const __m128d magF = _mm_load_sd(aligned_mags + f);

                            // Do the math
                            __m128d deltaMags = _mm_sub_sd(magF, obsF);   // magsF - obsF
                            deltaMags = _mm_mul_sd(deltaMags, deltaMags); // sqr (deltaMags)
                            deltaMags = _mm_mul_sd(deltaMags, varF);      // deltaMags * varF
                            deltaMags = _mm_add_sd(deltaMags, va2F);      // deltaMags + va2F
                            deltaMags = _mm_mul_sd(deltaMags, half);      // deltaMags * (-0.5)

                            likelihood = _mm_add_sd(likelihood, deltaMags);
                        }

                        const __m128d posts = Sleef_expd2_u10sse2(likelihood);

                        __m128d tPost = _mm_load_sd(sPost);

                        tPost = _mm_add_sd(tPost, posts);
                        _mm_store_sd(sPost, tPost);
                    }
                }
            }
        }
    }

    vector<double> vPost;
    vPost.assign(post, post + nSystems);

    _mm_free(aligned_mags);
    _mm_free(post);

    return vPost;
}

#elif BUILD_PARALLEL

// A (parallel) reference implementation of the above. This
// implementation does not require pre-allocated/aligned memory or any
// pre-computation, It also doesn't use SIMD and runs some 50x slower
// than the above. YMMV.

vector<double> margEvolveNoBinaries(const Cluster &clust, const Model &evoModels, const Isochrone &isochrone, ThreadPool &pool, vector<StellarSystem> &systems, const bool modIsParallax, const int eepInterpolationPower)
{
    mutex logPostMutex;
    const int isoIncrem = GlobalSettings::calcEepSteps(eepInterpolationPower);

    vector<double> secondaryMags;

    {
        Star s;
        s.mass = 0.0;

        secondaryMags = s.getMags(clust, evoModels, isochrone);
    }

    vector<double> post(systems.size(), 0.0);

    {
        auto agbTipMass = isochrone.agbTipMass();
        auto isoSize    = isochrone.eeps.size();
        auto nSystems   = systems.size();

        pool.parallelFor(isoSize - 1, [=,&post, &isochrone, &evoModels, &clust, &secondaryMags, &logPostMutex](int m)
        {
            double dIsoMass = isochrone.eeps[m + 1].mass - isochrone.eeps[m].mass;

            // In the event that we have an invalid range, skip that range
            // This generally occurs only at very high EEPs, where the masses are close together
            if (dIsoMass > 0.0)
            {
                Star s;
                vector<double> primaryMags, combinedMags, threadPost(nSystems, 0.0);

                double dMass    = dIsoMass / isoIncrem;
                double logdMass = __builtin_log (dMass);

                double cMass = isochrone.eeps[m].mass;

                for (auto k = 0; k < isoIncrem; ++k)
                {
                    double primaryMass = cMass + (k * dMass);

                    if (primaryMass >= agbTipMass)
                    {
                        primaryMass = agbTipMass;

                        // Kludgey way of exiting the for loop
                        // Only works if we don't use K below here
                        // Ensures we add the agbTip likelihood once and only once
                        k = isoIncrem;
                    }

                    const double logPrior    = clust.logPriorMass (primaryMass);

                    s.mass       = primaryMass;
                    primaryMags  = s.getMags(clust, evoModels, isochrone);
                    combinedMags = StellarSystem::deriveCombinedMags(clust, evoModels, isochrone, primaryMags, secondaryMags, modIsParallax);

                    for (size_t s = 0; s < nSystems; ++s)
                    {
                        threadPost[s] += __builtin_exp(systems[s].logPost (clust, evoModels, isochrone, logPrior, combinedMags) + logdMass);
                    }
                }

                std::lock_guard<mutex> lk(logPostMutex);

                for (size_t s = 0; s < nSystems; ++s)
                    post[s] += threadPost[s];
            }
        });
    }

    return post;
}

#else

#include <arm_neon.h>
#include <sleef.h>

#define ROUND_DOWN(x, s) ((x) & ~((s)-1))

inline static float64x2_t calcDeltaMags_neon(const double* sVar, const double* sVa2, const double* sObs, const float64x2_t magF, const float64x2_t half, const size_t f)
{
    const float64x2_t varF = vld1q_f64(sVar + f);
    const float64x2_t va2F = vld1q_f64(sVa2 + f);
    const float64x2_t obsF = vld1q_f64(sObs + f);

    // Do the math
    float64x2_t deltaMags = vsubq_f64(magF, obsF);    // magsF - obsF
    deltaMags = vmulq_f64(deltaMags, deltaMags);      // sqr (deltaMags)
    deltaMags = vmulq_f64(deltaMags, varF);           // deltaMags * varF
    deltaMags = vaddq_f64(deltaMags, va2F);           // deltaMags + va2F
    deltaMags = vmulq_f64(deltaMags, half);           // deltaMags * (-0.5)

    return deltaMags;
}

inline static float64x2_t calcSingleFilterResult_neon(const double* sVar1, const double* sVa21, const double* sObs1, const double* sVar2, const double* sVa22, const double* sObs2, const float64x2_t magF, const float64x2_t half, const size_t f)
{
    // Load low element from star 1, high from star 2
    float64x2_t varF = vld1q_lane_f64(sVar1 + f, vdupq_n_f64(0), 0);
    varF = vld1q_lane_f64(sVar2 + f, varF, 1);

    float64x2_t va2F = vld1q_lane_f64(sVa21 + f, vdupq_n_f64(0), 0);
    va2F = vld1q_lane_f64(sVa22 + f, va2F, 1);

    float64x2_t obsF = vld1q_lane_f64(sObs1 + f, vdupq_n_f64(0), 0);
    obsF = vld1q_lane_f64(sObs2 + f, obsF, 1);

    // Do the normal math for one filter across two stars
    float64x2_t deltaMags = vsubq_f64(magF, obsF);    // magsF - obsF
    deltaMags = vmulq_f64(deltaMags, deltaMags);      // sqr (deltaMags)
    deltaMags = vmulq_f64(deltaMags, varF);           // deltaMags * varF
    deltaMags = vaddq_f64(deltaMags, va2F);           // deltaMags + va2F
    deltaMags = vmulq_f64(deltaMags, half);           // deltaMags * (-0.5)

    return deltaMags;
}

inline static void addLogLikelihood_neon(const float64x2_t likelihood, double* sPost)
{
    const float64x2_t posts = Sleef_expd2_u10advsimd(likelihood);

    float64x2_t tPost = vld1q_f64(sPost);
    tPost = vaddq_f64(tPost, posts);
    vst1q_f64(sPost, tPost);
}

// ARM64 NEON optimized version
vector<double> margEvolveNoBinaries(const Cluster &clust, const Model &evoModels, const Isochrone &isochrone, ThreadPool &, const double* const systemVars, const double* const systemVar2, const double* const systemObs, const size_t nSystems, const size_t obsSize, const size_t obsUnaligned, const bool modIsParallax, const int eepInterpolationPower)
{
    const int isoIncrem = GlobalSettings::calcEepSteps(eepInterpolationPower);

    double* post = reinterpret_cast<double*>(aligned_alloc(16, nSystems * sizeof(double)));

    // Correct distance if it's in parallax
    double distance;

    if (modIsParallax)
    {
        auto parsecs = 1 / clust.mod;
        distance = clust.abs + 5 * log10(parsecs) - 5;
    }
    else
    {
        distance = clust.mod;
    }

    {
        size_t i;
        auto tPost = post;
        const float64x2_t zero = vdupq_n_f64(0.0);

        for (i = 0; i < ROUND_DOWN(nSystems, 8); i += 8, tPost += 8)
        {
            vst1q_f64(tPost,     zero);
            vst1q_f64(tPost + 2, zero);
            vst1q_f64(tPost + 4, zero);
            vst1q_f64(tPost + 6, zero);
        }
        for (; i < nSystems; ++i)
            post[i] = 0;
    }

    double* aligned_mags = reinterpret_cast<double*>(aligned_alloc(16, obsSize * sizeof(double)));

    {
        auto agbTipMass = isochrone.agbTipMass();
        auto isoSize    = isochrone.eeps.size();

        for (size_t m = 0; m < isoSize - 1; ++m)
        {
            double dIsoMass = isochrone.eeps[m + 1].mass - isochrone.eeps[m].mass;

            // In the event that we have an invalid range, skip that range
            // This generally occurs only at very high EEPs, where the masses are close together
            if (dIsoMass > 0.0)
            {
                Star s;
                vector<double> combinedMags;

                double dMass    = dIsoMass / isoIncrem;
                double logdMass = __builtin_log (dMass);

                double cMass = isochrone.eeps[m].mass;

                for (auto k = 0; k < isoIncrem; ++k)
                {
                    double primaryMass = cMass + (k * dMass);

                    if (primaryMass >= agbTipMass)
                    {
                        primaryMass = agbTipMass;

                        // Kludgey way of exiting the for loop
                        // Only works if we don't use K below here
                        // Ensures we don't add the agbTip likelihood more than once
                        k = isoIncrem;
                    }

                    const double logPrior = clust.logPriorMass (primaryMass);

                    s.mass       = primaryMass;
                    combinedMags = s.msRgbEvol(isochrone);

                    for (size_t f = 0; f < obsUnaligned; ++f)
                    {
                        combinedMags[f] += distance;
                        combinedMags[f] += (evoModels.absCoeffs[f] - 1.0) * clust.abs;
                    }

                    std::copy(combinedMags.begin(), combinedMags.end(), aligned_mags);

                    // Pointers to first 8 stars
                    auto sVar1 = systemVars;
                    auto sVa21 = systemVar2;
                    auto sObs1 = systemObs;

                    auto sVar2 = sVar1 + obsSize;
                    auto sVa22 = sVa21 + obsSize;
                    auto sObs2 = sObs1 + obsSize;

                    auto sVar3 = sVar2 + obsSize;
                    auto sVa23 = sVa22 + obsSize;
                    auto sObs3 = sObs2 + obsSize;

                    auto sVar4 = sVar3 + obsSize;
                    auto sVa24 = sVa23 + obsSize;
                    auto sObs4 = sObs3 + obsSize;

                    auto sVar5 = sVar4 + obsSize;
                    auto sVa25 = sVa24 + obsSize;
                    auto sObs5 = sObs4 + obsSize;

                    auto sVar6 = sVar5 + obsSize;
                    auto sVa26 = sVa25 + obsSize;
                    auto sObs6 = sObs5 + obsSize;

                    auto sVar7 = sVar6 + obsSize;
                    auto sVa27 = sVa26 + obsSize;
                    auto sObs7 = sObs6 + obsSize;

                    auto sVar8 = sVar7 + obsSize;
                    auto sVa28 = sVa27 + obsSize;
                    auto sObs8 = sObs7 + obsSize;

                    auto sPost = post;

                    size_t st;

                    // Process 8 stars at a time (4 pairs of 2)
                    const auto perIter = 8;
                    const auto skip = obsSize * perIter;

                    const float64x2_t half   = vdupq_n_f64(-0.5);
                    const float64x2_t priors = vdupq_n_f64(logPrior + logdMass);

                    for (st = 0; st < ROUND_DOWN(nSystems, perIter); st += perIter, sPost += perIter,
                             sVar1 += skip, sVa21 += skip, sObs1 += skip,
                             sVar2 += skip, sVa22 += skip, sObs2 += skip,
                             sVar3 += skip, sVa23 += skip, sObs3 += skip,
                             sVar4 += skip, sVa24 += skip, sObs4 += skip,
                             sVar5 += skip, sVa25 += skip, sObs5 += skip,
                             sVar6 += skip, sVa26 += skip, sObs6 += skip,
                             sVar7 += skip, sVa27 += skip, sObs7 += skip,
                             sVar8 += skip, sVa28 += skip, sObs8 += skip)
                    {
                        float64x2_t likelihood1 = priors;
                        float64x2_t likelihood2 = priors;
                        float64x2_t likelihood3 = priors;
                        float64x2_t likelihood4 = priors;

                        size_t f;

                        // Main loop for filters: 8 stars, 2 filters per iteration
                        for (f = 0; f < ROUND_DOWN(obsUnaligned, 2); f += 2)
                        {
                            // Shared magnitudes for this filter pair
                            const float64x2_t magF = vld1q_f64(aligned_mags + f);

                            // Calculate delta magnitudes for all 8 stars
                            float64x2_t deltaMags1 = calcDeltaMags_neon(sVar1, sVa21, sObs1, magF, half, f);
                            float64x2_t deltaMags2 = calcDeltaMags_neon(sVar2, sVa22, sObs2, magF, half, f);
                            float64x2_t deltaMags3 = calcDeltaMags_neon(sVar3, sVa23, sObs3, magF, half, f);
                            float64x2_t deltaMags4 = calcDeltaMags_neon(sVar4, sVa24, sObs4, magF, half, f);
                            float64x2_t deltaMags5 = calcDeltaMags_neon(sVar5, sVa25, sObs5, magF, half, f);
                            float64x2_t deltaMags6 = calcDeltaMags_neon(sVar6, sVa26, sObs6, magF, half, f);
                            float64x2_t deltaMags7 = calcDeltaMags_neon(sVar7, sVa27, sObs7, magF, half, f);
                            float64x2_t deltaMags8 = calcDeltaMags_neon(sVar8, sVa28, sObs8, magF, half, f);

                            // Horizontal add pairs of filter results
                            const float64x2_t result1 = vpaddq_f64(deltaMags1, deltaMags2);
                            const float64x2_t result2 = vpaddq_f64(deltaMags3, deltaMags4);
                            const float64x2_t result3 = vpaddq_f64(deltaMags5, deltaMags6);
                            const float64x2_t result4 = vpaddq_f64(deltaMags7, deltaMags8);

                            // Accumulate likelihoods
                            likelihood1 = vaddq_f64(likelihood1, result1);
                            likelihood2 = vaddq_f64(likelihood2, result2);
                            likelihood3 = vaddq_f64(likelihood3, result3);
                            likelihood4 = vaddq_f64(likelihood4, result4);
                        }

                        // Handle odd number of filters
                        for (; f < obsUnaligned; ++f)
                        {
                            // Broadcast single magnitude
                            const float64x2_t magF = vdupq_n_f64(aligned_mags[f]);

                            const float64x2_t result1 = calcSingleFilterResult_neon(sVar1, sVa21, sObs1, sVar2, sVa22, sObs2, magF, half, f);
                            const float64x2_t result2 = calcSingleFilterResult_neon(sVar3, sVa23, sObs3, sVar4, sVa24, sObs4, magF, half, f);
                            const float64x2_t result3 = calcSingleFilterResult_neon(sVar5, sVa25, sObs5, sVar6, sVa26, sObs6, magF, half, f);
                            const float64x2_t result4 = calcSingleFilterResult_neon(sVar7, sVa27, sObs7, sVar8, sVa28, sObs8, magF, half, f);

                            likelihood1 = vaddq_f64(likelihood1, result1);
                            likelihood2 = vaddq_f64(likelihood2, result2);
                            likelihood3 = vaddq_f64(likelihood3, result3);
                            likelihood4 = vaddq_f64(likelihood4, result4);
                        }

                        addLogLikelihood_neon(likelihood1, sPost);
                        addLogLikelihood_neon(likelihood2, sPost + 2);
                        addLogLikelihood_neon(likelihood3, sPost + 4);
                        addLogLikelihood_neon(likelihood4, sPost + 6);
                    }

                    // Cleanup loop for remaining stars (< 8)
                    for (; st < nSystems; ++st, ++sPost,
                             sVar1 += obsSize, sVa21 += obsSize, sObs1 += obsSize)
                    {
                        double likelihood = logPrior + logdMass;

                        for (size_t f = 0; f < obsUnaligned; ++f)
                        {
                            double varF = sVar1[f];
                            double va2F = sVa21[f];
                            double obsF = sObs1[f];
                            double magF = aligned_mags[f];

                            double deltaMags = magF - obsF;
                            deltaMags = deltaMags * deltaMags;
                            deltaMags = deltaMags * varF;
                            deltaMags = deltaMags + va2F;
                            deltaMags = deltaMags * (-0.5);

                            likelihood += deltaMags;
                        }

                        *sPost += __builtin_exp(likelihood);
                    }
                }
            }
        }
    }

    vector<double> vPost;
    vPost.assign(post, post + nSystems);

    free(aligned_mags);
    free(post);

    return vPost;
}

#endif

// Scalar reference implementation - always compiled for testing equivalence.
// This is the canonical implementation against which optimized versions are verified.
vector<double> margEvolveNoBinaries_Reference(const Cluster &clust, const Model &evoModels, const Isochrone &isochrone, vector<StellarSystem> &systems, const bool modIsParallax, const int eepInterpolationPower)
{
    const int isoIncrem = GlobalSettings::calcEepSteps(eepInterpolationPower);

    vector<double> secondaryMags;

    {
        Star s;
        s.mass = 0.0;

        secondaryMags = s.getMags(clust, evoModels, isochrone);
    }

    vector<double> post(systems.size(), 0.0);

    {
        auto agbTipMass = isochrone.agbTipMass();
        auto isoSize    = isochrone.eeps.size();
        auto nSystems   = systems.size();

        for (size_t m = 0; m < isoSize - 1; ++m)
        {
            double dIsoMass = isochrone.eeps[m + 1].mass - isochrone.eeps[m].mass;

            // In the event that we have an invalid range, skip that range
            // This generally occurs only at very high EEPs, where the masses are close together
            if (dIsoMass > 0.0)
            {
                Star s;
                vector<double> primaryMags, combinedMags;

                double dMass    = dIsoMass / isoIncrem;
                double logdMass = __builtin_log (dMass);

                double cMass = isochrone.eeps[m].mass;

                for (auto k = 0; k < isoIncrem; ++k)
                {
                    double primaryMass = cMass + (k * dMass);

                    if (primaryMass >= agbTipMass)
                    {
                        primaryMass = agbTipMass;

                        // Kludgey way of exiting the for loop
                        // Only works if we don't use K below here
                        // Ensures we add the agbTip likelihood once and only once
                        k = isoIncrem;
                    }

                    const double logPrior = clust.logPriorMass (primaryMass);

                    s.mass       = primaryMass;
                    primaryMags  = s.getMags(clust, evoModels, isochrone);
                    combinedMags = StellarSystem::deriveCombinedMags(clust, evoModels, isochrone, primaryMags, secondaryMags, modIsParallax);

                    for (size_t s = 0; s < nSystems; ++s)
                    {
                        post[s] += __builtin_exp(systems[s].logPost (clust, evoModels, isochrone, logPrior, combinedMags) + logdMass);
                    }
                }
            }
        }
    }

    return post;
}
