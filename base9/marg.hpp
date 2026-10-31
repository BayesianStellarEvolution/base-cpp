#ifndef MARG_HPP
#define MARG_HPP

#include <vector>

#include "Star.hpp"
#include "Cluster.hpp"

#include "Model.hpp"
#include "Utility.hpp"

std::vector<double> margEvolveWithBinary (const Cluster &, std::vector<StellarSystem>&, const Model&, const Isochrone&, base::utility::ThreadPool&, const bool, const int);

// SIMD-optimized binary version - takes pre-aligned arrays
// Additional arrays: systemFluxLo (exp10((obs - 3*sqrt(var)) / -2.5)), systemFluxHi (exp10((obs + 3*sqrt(var)) / -2.5))
std::vector<double> margEvolveWithBinary(const Cluster&, const Model&, const Isochrone&, base::utility::ThreadPool&, std::vector<StellarSystem>&, const double* const, const double* const, const double* const, const double* const, const double* const, const size_t, const size_t, const bool, const int);

// SIMD-optimized version (x86 SSE or ARM NEON) - takes pre-aligned arrays
std::vector<double> margEvolveNoBinaries(const Cluster&, const Model&, const Isochrone&, base::utility::ThreadPool&, const double* const, const double* const, const double* const, const size_t, const size_t, const size_t, const bool, const int);

// Scalar reference implementation - always available for testing
std::vector<double> margEvolveNoBinaries_Reference(const Cluster &, const Model &, const Isochrone &, std::vector<StellarSystem> &, const bool, const int);

#endif
