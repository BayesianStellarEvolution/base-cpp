#include <array>
#include <atomic>
#include <mutex>

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "Cluster.hpp"
#include "Star.hpp"

#include "constants.hpp"
#include "evolve.hpp"
#include "structures.hpp"
#include "densities.hpp"

using std::array;
using std::atomic;
using std::once_flag;
using std::call_once;

constexpr double sqr(double a)
{
    return a * a;
}

static_assert(M_PI > 3.141592, "M_PI is defined and and at least 3.141592");
static_assert(M_PI < 3.15, "M_PI is defined and less than 3.15");

const double mf_sigma = 0.67729, mf_mu = -1.02;
const double loglog10 = log (log (10));

namespace evil
{
    // Calculate the mass normalization factor once so that we don't have to
    // calculate it every time we want the mass prior
    class logMassNorm
    {
      private:
        logMassNorm(double M_wd_up)
        {
            double p, q, c;
            double tup, tlow;

            p = mf_mu + mf_sigma * mf_sigma * log (10);
            tup = (log10 (M_wd_up) - p) / (mf_sigma);
            tlow = (-1 - p) / mf_sigma;
            q = exp (-(mf_mu * mf_mu - p * p) / (2 * mf_sigma * mf_sigma));
            c = 1 / (q * mf_sigma * sqrt (2 * M_PI) * (Phi (tup) - Phi (tlow)));

            this->var = log (c);
        }
        logMassNorm(const logMassNorm&) = delete;
        logMassNorm(const logMassNorm&&) = delete;
        logMassNorm& operator=(const logMassNorm&) = delete;

        double var;

      public:
        static logMassNorm& getInstance(const Cluster &pCluster)
        {
            static logMassNorm instance(pCluster.M_wd_up);

            return instance;
        }

        double getLogMassNorm() const { return var; }
    };
}

double logPriorMass (const Star &pStar, const Cluster &pCluster)
// Compute log prior density
{
    double mass1, log_m1, logPrior = 0.0;

    mass1 = pStar.getMass1();

    assert (mass1 > 0.1 && mass1 <= pCluster.M_wd_up);

    if (pStar.isFieldStar)
    {
        logPrior = logTDens (pStar.U, pStar.meanU, pStar.varU, DOF);
        if (pStar.status.at(0) != WD)
            logPrior += logTDens (pStar.massRatio, pStar.meanMassRatio, pStar.varMassRatio, DOF);
        return logPrior;
    }
    else
    {
        log_m1 = log10 (mass1);
        logPrior = evil::logMassNorm::getInstance(pCluster).getLogMassNorm() + -0.5 * sqr (log_m1 - mf_mu) / (sqr (mf_sigma)) - log (mass1) - loglog10;
        return logPrior;
    }
}

// Compute log prior density for cluster properties
double logPriorClust (const Cluster &pCluster, const Model &evoModels)
{
    if ((pCluster.age < evoModels.mainSequenceEvol->getMinAge())
        || (pCluster.age > evoModels.mainSequenceEvol->getMaxAge())
        || (pCluster.ifmrSlope < 0.0)
        || (pCluster.abs < 0.0)
        || ((evoModels.IFMR == 11) && (pCluster.ifmrQuadCoef < 0.0)))
    {
        throw InvalidCluster("Invalid cluster parameter");
    }

    // enforce monotonicity in IFMR
    if (evoModels.IFMR == 10)
    {
        double massLower = 0.15;
        double massUpper = pCluster.M_wd_up;
        double massShift = 3.0;
        double angle = atan (pCluster.ifmrSlope);
        double aa = cos (angle) * (1 + pCluster.ifmrSlope * pCluster.ifmrSlope);
        double xLower = aa * (massLower - massShift);
        double xUpper = aa * (massUpper - massShift);

        double dydx_xLower = pCluster.ifmrQuadCoef * (xLower - xUpper);
        double dydx_xUpper = -dydx_xLower;

        double slopeLower = tan (angle + atan (dydx_xLower));
        double slopeUpper = tan (angle + atan (dydx_xUpper));

        // if IFMR is decreasing at either endpoint, reject
        if (slopeLower < 0.0 || slopeUpper < 0.0)
            throw std::range_error("IFMR decreasing at endpoint");
    }

    double prior = 0.0;
    //DS: with a uniform prior on carbonicity, the above won't change since log(1) = 0.

    if (pCluster.priorVar.at(FEH) > EPSILON)
        prior += (-0.5) * sqr (pCluster.feh - pCluster.priorMean.at(FEH)) / pCluster.priorVar.at(FEH);
    if (pCluster.priorVar.at(MOD) > EPSILON)
        prior += (-0.5) * sqr (pCluster.mod - pCluster.priorMean.at(MOD)) / pCluster.priorVar.at(MOD);
    if (pCluster.priorVar.at(ABS) > EPSILON)
        prior += (-0.5) * sqr (pCluster.abs - pCluster.priorMean.at(ABS)) / pCluster.priorVar.at(ABS);
    if (pCluster.priorVar.at(YYY) > EPSILON)
        prior += (-0.5) * sqr (pCluster.yyy - pCluster.priorMean.at(YYY)) / pCluster.priorVar.at(YYY);

    return prior;
}

double logLikelihood (int numFilts, const Star &pStar, const array<double, FILTS> &mags, const array<double, FILTS> &filterPriorMin, const array<double, FILTS> &filterPriorMax)
// Computes log likelihood
{
    int i;
    double likelihood = 0.0;

    for (i = 0; i < numFilts; i++)
    {
        if (pStar.isFieldStar)
        {
            assert (filterPriorMin.at(i) <= pStar.obsPhot.at(i) && pStar.obsPhot.at(i) <= filterPriorMax.at(i));

            likelihood -= log (filterPriorMax.at(i) - filterPriorMin.at(i));
        }
        else
        {
            if (pStar.variance.at(i) > 1e-9)
                likelihood -= 0.5 * (log (2 * M_PI * pStar.variance.at(i)) + (sqr (mags.at(i) - pStar.obsPhot.at(i)) / pStar.variance.at(i)));
        }
    }
    return likelihood;
}

double tLogLikelihood (int numFilts, const Star &pStar, const array<double, FILTS> &mags, const array<double, FILTS> &filterPriorMin, const array<double, FILTS> &filterPriorMax)
// Computes log likelihood
{
    int i;
    double likelihood = 0.0;
    double dof = 3.0;
    double quadratic_sum = 0.0;

    if (pStar.isFieldStar)
    {
        for (i = 0; i < numFilts; i++)
        {
            assert (filterPriorMin.at(i) <= pStar.obsPhot.at(i) && pStar.obsPhot.at(i) <= filterPriorMax.at(i));

            likelihood -= log (filterPriorMax.at(i) - filterPriorMin.at(i));
        }
    }
    else
    {
        for (i = 0; i < numFilts; i++)
        {
            if (pStar.variance.at(i) > 1e-9)
            {
                quadratic_sum += sqr (mags.at(i) - pStar.obsPhot.at(i)) / pStar.variance.at(i);
                likelihood -= 0.5 * (log (M_PI * pStar.variance.at(i)));
            }
        }
        likelihood += lgamma (0.5 * (dof + (double) numFilts)) - lgamma (0.5 * dof);
        likelihood -= 0.5 * (dof + (double) numFilts) * log (1 + quadratic_sum);
    }
    return likelihood;
}

double scaledLogLike (int numFilts, const Star &pStar, double varScale, const array<double, FILTS> &mags, const array<double, FILTS> &filterPriorMin, const array<double, FILTS> &filterPriorMax)
// Computes log likelihood
{
    int i;
    double likelihood = 0.0;

    for (i = 0; i < numFilts; i++)
    {
        if (pStar.isFieldStar)
        {
            assert (filterPriorMin.at(i) <= pStar.obsPhot.at(i) && pStar.obsPhot.at(i) <= filterPriorMax.at(i));

            likelihood -= log (filterPriorMax.at(i) - filterPriorMin.at(i));
        }
        else
        {
            if (pStar.variance.at(i) > 1e-9)
            {
                likelihood -= 0.5 * (log (2 * M_PI * varScale * pStar.variance.at(i)) + (sqr (mags.at(i) - pStar.obsPhot.at(i)) / (varScale * pStar.variance.at(i))));
            }

        }
    }
    return likelihood;
}


double logPost1Star (const Star &pStar, const Cluster &pCluster, const Model &evoModels, const array<double, FILTS> &mags, const array<double, FILTS> &filterPriorMin, const array<double, FILTS> &filterPriorMax)
// Compute posterior density for 1 star:
{
    double likelihood = 0.0, logPrior = 0.0;

    logPrior = logPriorMass (pStar, pCluster);

    likelihood = scaledLogLike (evoModels.numFilts, pStar, pCluster.varScale, mags, filterPriorMin, filterPriorMax);

    return (logPrior + likelihood);
}


// computes normal distribution Phi(x) (integral from -Inf to x of normal density)
// taken from: http://www.jstatsoft.org/v11/i04/v11i04.pdf
double Phi (double x)
{
    long double s = x, t = 0, b = x, q = x * x, i = 1;

    while (s != t)
        s = (t = s) + (b *= q / (i += 2));

    return 0.5 + s * exp (-0.5 * q - 0.91893853320467274178L);
}

// do not call this routine with nu = 2, which wouldn't make much sense anyway
double logTDens (double x, double mean, double var, double nu)
{
    double logp = 0;
    double s;

    s = sqrt (nu / (var * (nu - 2)));

    logp = log (s) + gamma(nu) - 3.5 * log (1 + pow (s * (x - mean), 2) / nu);

    return logp;
}
