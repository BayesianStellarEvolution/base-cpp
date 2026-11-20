// Copyright (c) Embry-Riddle Aeronautical University. All rights reserved.
//
// Compare log-likelihood values across different EEP interpolation powers
// to assess whether lower powers preserve the likelihood surface shape.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
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

struct LikelihoodResult {
    double logLike;
    int skippedStars;
};

LikelihoodResult computeLogLikelihood(
    Cluster& clust,
    Model& evoModels,
    Isochrone& isochrone,
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
    vector<double> post;
    if (noBinaries) {
        post = margEvolveNoBinaries(clust, evoModels, isochrone, pool,
                                    sysVars, sysVar2, sysObs,
                                    systems.size(), howManyFiltsAligned, howManyFilts,
                                    false, eepPower);
    } else {
        post = margEvolveWithBinary(clust, evoModels, isochrone, pool,
                                    systems,
                                    sysVars, sysVar2, sysObs, sysFluxLo, sysFluxHi,
                                    howManyFiltsAligned, howManyFilts,
                                    false, eepPower);
    }

    double logLike = 0.0;
    int skipped = 0;
    for (size_t i = 0; i < systems.size(); ++i) {
        if (post[i] > 0.0) {
            logLike += std::log(post[i]);
        } else {
            ++skipped;
        }
    }

    return {logLike, skipped};
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
    Model evoModels = makeModel(settings.global);

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
    const double ageCenter = baseClust.age;
    const double ageSpan = 0.2;
    const double fehCenter = baseClust.feh;
    const double fehSpan = 0.2;

    const int nPowers = settings.maxPower - settings.minPower + 1;
    const int nGridPoints = settings.gridPoints * settings.gridPoints;
    const bool refInRange = (settings.referencePower >= settings.minPower &&
                             settings.referencePower <= settings.maxPower);

    // Pre-compute grid coordinates (age x feh, distMod fixed - it doesn't interact with EEP power)
    const double gridDenom = 1.0 / (settings.gridPoints - 1);
    vector<double> ageGrid(settings.gridPoints), fehGrid(settings.gridPoints);
    for (int i = 0; i < settings.gridPoints; ++i) {
        double frac = i * gridDenom;
        ageGrid[i] = ageCenter - ageSpan * 0.5 + ageSpan * frac;
        fehGrid[i] = fehCenter - fehSpan * 0.5 + fehSpan * frac;
    }

    vector<vector<double>> logLikes(nPowers, vector<double>(nGridPoints));
    vector<double> refLogLikes(nGridPoints);
    vector<long> skippedStars(nPowers, 0);
    long refSkippedStars = 0;

    cerr << "Computing likelihoods on " << nGridPoints << " grid points for powers "
         << settings.minPower << " to " << settings.maxPower;
    if (!refInRange) {
        cerr << " (plus reference power " << settings.referencePower << ")";
    }
    cerr << "..." << endl;

    // Reusable cluster for likelihood computation (marg functions don't modify it)
    Cluster clust = baseClust;

    int gridIdx = 0;
    for (int iAge = 0; iAge < settings.gridPoints; ++iAge) {
        const double age = ageGrid[iAge];
        clust.age = age;

        for (int iFeh = 0; iFeh < settings.gridPoints; ++iFeh) {
            const double feh = fehGrid[iFeh];
            clust.feh = feh;

            // Derive isochrone once per (age, feh) pair
            Isochrone* isochrone = evoModels.mainSequenceEvol->deriveIsochrone(feh, clust.yyy, age);

            for (int power = settings.minPower; power <= settings.maxPower; ++power) {
                auto result = computeLogLikelihood(
                    clust, evoModels, *isochrone, systems,
                    sysVars.data(), sysVar2.data(), sysObs.data(),
                    sysFluxLo.data(), sysFluxHi.data(),
                    howManyFiltsAligned, howManyFilts,
                    pool, settings.noBinaries, power);
                logLikes[power - settings.minPower][gridIdx] = result.logLike;
                skippedStars[power - settings.minPower] += result.skippedStars;
            }

            if (refInRange) {
                refLogLikes[gridIdx] = logLikes[settings.referencePower - settings.minPower][gridIdx];
            } else {
                auto result = computeLogLikelihood(
                    clust, evoModels, *isochrone, systems,
                    sysVars.data(), sysVar2.data(), sysObs.data(),
                    sysFluxLo.data(), sysFluxHi.data(),
                    howManyFiltsAligned, howManyFilts,
                    pool, settings.noBinaries, settings.referencePower);
                refLogLikes[gridIdx] = result.logLike;
                refSkippedStars += result.skippedStars;
            }
            ++gridIdx;

            delete isochrone;
        }
        cerr << "  Age slice " << (iAge + 1) << "/" << settings.gridPoints << " complete" << endl;
    }

    cout << "\n=== EEP Power Comparison Results ===" << endl;
    cout << "Reference power: " << settings.referencePower << endl;
    cout << "Grid points: " << nGridPoints << " (" << settings.gridPoints << "^2 age x feh)" << endl;
    cout << "Grid spacing: age=" << std::fixed << std::setprecision(4) << (ageSpan / (settings.gridPoints - 1))
         << ", feh=" << (fehSpan / (settings.gridPoints - 1)) << endl;
    cout << "Stars: " << systems.size() << ", Filters: " << howManyFilts << endl;

    // Show mass spacing for each power using isochrone at center values
    cout << "\n--- Mass Resolution per EEP Power ---" << endl;
    {
        Isochrone* sampleIso = evoModels.mainSequenceEvol->deriveIsochrone(fehCenter, baseClust.yyy, ageCenter);
        cout << std::setw(8) << "Power"
             << std::setw(10) << "Steps"
             << std::setw(14) << "Min dMass"
             << std::setw(14) << "Max dMass"
             << std::setw(14) << "Mean dMass"
             << endl;
        cout << std::string(60, '-') << endl;

        for (int power = settings.minPower; power <= settings.maxPower; ++power) {
            int isoIncrem = GlobalSettings::calcEepSteps(power);
            double minDMass = std::numeric_limits<double>::max();
            double maxDMass = 0.0;
            double sumDMass = 0.0;
            int nSegments = 0;

            for (size_t m = 0; m < sampleIso->eeps.size() - 1; ++m) {
                double dIsoMass = sampleIso->eeps[m + 1].mass - sampleIso->eeps[m].mass;
                if (dIsoMass > 0.0) {
                    double dMass = dIsoMass / isoIncrem;
                    if (dMass < minDMass) minDMass = dMass;
                    if (dMass > maxDMass) maxDMass = dMass;
                    sumDMass += dMass;
                    ++nSegments;
                }
            }

            cout << std::setw(8) << power
                 << std::setw(10) << isoIncrem
                 << std::setw(14) << std::scientific << std::setprecision(4) << minDMass
                 << std::setw(14) << maxDMass
                 << std::setw(14) << (nSegments > 0 ? sumDMass / nSegments : 0.0)
                 << endl;
        }

        // Include reference power if out of range
        if (!refInRange) {
            int isoIncrem = GlobalSettings::calcEepSteps(settings.referencePower);
            double minDMass = std::numeric_limits<double>::max();
            double maxDMass = 0.0;
            double sumDMass = 0.0;
            int nSegments = 0;

            for (size_t m = 0; m < sampleIso->eeps.size() - 1; ++m) {
                double dIsoMass = sampleIso->eeps[m + 1].mass - sampleIso->eeps[m].mass;
                if (dIsoMass > 0.0) {
                    double dMass = dIsoMass / isoIncrem;
                    if (dMass < minDMass) minDMass = dMass;
                    if (dMass > maxDMass) maxDMass = dMass;
                    sumDMass += dMass;
                    ++nSegments;
                }
            }

            cout << std::setw(8) << settings.referencePower << " (ref)"
                 << std::setw(10) << isoIncrem
                 << std::setw(14) << std::scientific << std::setprecision(4) << minDMass
                 << std::setw(14) << maxDMass
                 << std::setw(14) << (nSegments > 0 ? sumDMass / nSegments : 0.0)
                 << endl;
        }

        delete sampleIso;
    }

    // Diagnostic: show raw log-likelihood values for each power
    cout << "\n--- Raw Log-Likelihood Summary ---" << endl;
    cout << std::setw(8) << "Power"
         << std::setw(14) << "Mean LL"
         << std::setw(14) << "Min LL"
         << std::setw(14) << "Max LL"
         << std::setw(12) << "Skipped"
         << endl;
    cout << std::string(62, '-') << endl;

    for (int power = settings.minPower; power <= settings.maxPower; ++power) {
        const int powerIdx = power - settings.minPower;
        const double* powerLikes = logLikes[powerIdx].data();

        double sum = 0.0;
        double minLL = std::numeric_limits<double>::max();
        double maxLL = std::numeric_limits<double>::lowest();

        for (int i = 0; i < nGridPoints; ++i) {
            sum += powerLikes[i];
            if (powerLikes[i] < minLL) minLL = powerLikes[i];
            if (powerLikes[i] > maxLL) maxLL = powerLikes[i];
        }

        cout << std::setw(8) << power
             << std::setw(14) << std::scientific << std::setprecision(4) << (sum / nGridPoints)
             << std::setw(14) << minLL
             << std::setw(14) << maxLL
             << std::setw(12) << std::fixed << std::setprecision(0) << skippedStars[powerIdx]
             << endl;
    }

    // Show reference separately if out of range
    if (!refInRange) {
        double sum = 0.0;
        double minLL = std::numeric_limits<double>::max();
        double maxLL = std::numeric_limits<double>::lowest();
        for (int i = 0; i < nGridPoints; ++i) {
            sum += refLogLikes[i];
            if (refLogLikes[i] < minLL) minLL = refLogLikes[i];
            if (refLogLikes[i] > maxLL) maxLL = refLogLikes[i];
        }
        cout << std::setw(8) << settings.referencePower << " (ref)"
             << std::setw(14) << std::scientific << std::setprecision(4) << (sum / nGridPoints)
             << std::setw(14) << minLL
             << std::setw(14) << maxLL
             << std::setw(12) << std::fixed << std::setprecision(0) << refSkippedStars
             << endl;
    }

    // Sample: show first few grid points for all powers
    cout << "\n--- Sample Grid Points (first 3) ---" << endl;
    cout << std::setw(6) << "Grid";
    for (int power = settings.minPower; power <= settings.maxPower; ++power) {
        cout << std::setw(14) << ("P" + std::to_string(power));
    }
    cout << endl;
    cout << std::string(6 + 14 * nPowers, '-') << endl;

    for (int i = 0; i < std::min(3, nGridPoints); ++i) {
        cout << std::setw(6) << i;
        for (int power = settings.minPower; power <= settings.maxPower; ++power) {
            cout << std::setw(14) << std::scientific << std::setprecision(4)
                 << logLikes[power - settings.minPower][i];
        }
        cout << endl;
    }

    cout << "\n--- Differences from Reference ---" << endl;
    cout << std::setw(8) << "Power"
         << std::setw(14) << "Mean Diff"
         << std::setw(14) << "Std Dev"
         << std::setw(14) << "Min Diff"
         << std::setw(14) << "Max Diff"
         << std::setw(10) << "Shape OK?"
         << endl;
    cout << std::string(74, '-') << endl;

    for (int power = settings.minPower; power <= settings.maxPower; ++power) {
        const int powerIdx = power - settings.minPower;
        const double* powerLikes = logLikes[powerIdx].data();
        const double* refLikes = refLogLikes.data();

        // Single-pass computation of sum, sumSq, min, max
        double sum = 0.0, sumSq = 0.0;
        double minDiff = std::numeric_limits<double>::max();
        double maxDiff = std::numeric_limits<double>::lowest();

        for (int i = 0; i < nGridPoints; ++i) {
            const double diff = powerLikes[i] - refLikes[i];
            sum += diff;
            sumSq += diff * diff;
            if (diff < minDiff) minDiff = diff;
            if (diff > maxDiff) maxDiff = diff;
        }

        const double mean = sum / nGridPoints;
        // Var(X) = E[X^2] - E[X]^2
        const double variance = sumSq / nGridPoints - mean * mean;
        const double stdDev = std::sqrt(std::max(0.0, variance));

        const bool shapeOK = (stdDev < 0.05);

        cout << std::setw(8) << power
             << std::setw(14) << std::scientific << std::setprecision(4) << mean
             << std::setw(14) << stdDev
             << std::setw(14) << minDiff
             << std::setw(14) << maxDiff
             << std::setw(10) << (shapeOK ? "YES" : "NO")
             << endl;
    }

    cout << "\nInterpretation:" << endl;
    cout << "  Mean Diff: constant offset (OK - cancels in MCMC ratios)" << endl;
    cout << "  Std Dev: shape distortion (BAD - changes posteriors)" << endl;
    cout << "  Shape OK: Std Dev < 0.05" << endl;

    return 0;
}
