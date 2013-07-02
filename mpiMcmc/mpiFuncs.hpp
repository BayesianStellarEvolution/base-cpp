#ifndef MPIFUNC_HPP
#define MPIFUNC_HPP

#include <array>
#include <fstream>

#include "Model.hpp"
#include "mpiMcmc.hpp"

void make_cholesky_decomp(struct ifmrMcmcControl &ctrl, Matrix<double, NPARAMS, nSave> &params);
double logPostStep(Chain &mc, const Model &evoModel, std::array<double, N_WD_MASS1> &wdMass1Grid, Cluster &propClust, double fsLike, std::array<double, 2> &ltau);
int acceptClustMarg (double logPostCurr, double logPostProp, std::array<double, 2> &ltau);
void initChain (Chain &mc, const struct ifmrMcmcControl &ctrl, const Model &evoModels, std::array<double, 2> &ltau);
void initIfmrMcmcControl (Chain &mc, struct ifmrMcmcControl &ctrl, const Model &evoModels, Settings &settings);
void initMassGrids (std::array<double, N_MS_MASS1 * N_MS_MASS_RATIO> &msMass1Grid, std::array<double, N_MS_MASS1 * N_MS_MASS_RATIO> &msMassRatioGrid, std::array<double, N_WD_MASS1> &wdMass1Grid, Chain const &mc);
void propClustBigSteps (Cluster &clust, struct ifmrMcmcControl const &ctrl);
void propClustIndep (Cluster &clust, struct ifmrMcmcControl const &ctrl);
void propClustCorrelated (Cluster &clust, struct ifmrMcmcControl const &ctrl);
void printHeader (std::ofstream &file, std::array<double, NPARAMS> const &priors);

#endif
