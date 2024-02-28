#include <cmath>
#include <iostream>
#include <random>
#include <gsl/gsl_sf_gamma.h>

#include "McmcApplication.hpp"

using std::cout;
using std::endl;

constexpr double sqr(double a)
{
    return a * a;
}

void beta_cdfs(int a, int b)
{
//  auto mid = acceptableCumulative(); // ?
    auto mid = gsl_sf_beta_inc(a, b, 0.4) - gsl_sf_beta_inc(a, b, 0.2);

    cout << gsl_sf_beta_inc(a, b, 0.2) << ", ";
//    std::cout << mid << ", " << (1 - mid) << std::endl;
    cout << mid << ", ";
    cout << 1 - gsl_sf_beta_inc(a, b, 0.4) << " - ";
    cout << sqrt((double)(a * b) / ((a + b) * (a + b) * (a + b + 1))) << endl; // Std. Dev
}

double McmcApplication::stdDev() const
{
    auto a = 1 + accepted;
    auto b = 1 + rejected;

    double top = a * b;
    double bot = sqr(a + b) * (a + b + 1);
    double res = top / bot;
    double fin = sqrt(res);

    if (std::isnan(fin))
        cout << top << ", " << bot << ", " << res << ", " << fin << endl;

    return fin;
}

double McmcApplication::lowerCumulative() const
{
    auto a = 1 + accepted;
    auto b = 1 + rejected;

    return gsl_sf_beta_inc(a, b, 0.2);
}

double McmcApplication::upperCumulative() const
{
    auto a = 1 + accepted;
    auto b = 1 + rejected;

    return 1 - gsl_sf_beta_inc(a, b, 0.4);
}

double McmcApplication::acceptableCumulative() const
{
    auto a = 1 + accepted;
    auto b = 1 + rejected;

    return gsl_sf_beta_inc(a, b, 0.4)
         - gsl_sf_beta_inc(a, b, 0.2);
}

void McmcApplication::reportAcceptanceRatio() const
{
    cout << acceptanceRatio() << ": ";
    cout << lowerCumulative() << ", ";
    cout << acceptableCumulative() << ", ";
    cout << upperCumulative() << " - ";
    cout << stdDev() << endl;
}

// Decides whether to accept a proposed cluster property
bool McmcApplication::acceptClustMarg (const double logPostCurr, const double logPostProp)
{
    if (std::isinf (logPostProp))
    {
//        cerr << "-Inf posterior proposed and rejected" << endl;
        rejected += 1;
        return false;
    }

    double alpha = logPostProp - logPostCurr;

    if (alpha >= 0)             // Short circuit exit to the MH algorithm
    {
        accepted += 1;
        return true;
    }

    double u = std::generate_canonical<double, 53>(gen);

    if (u < 1.e-15)
        u = 1.e-15;
    u = log (u);

    if (u < alpha)
    {
        accepted += 1;
        return true;
    }
    else
    {
        rejected += 1;
        return false;
    }
}
