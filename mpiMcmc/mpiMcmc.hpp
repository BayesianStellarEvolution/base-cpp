#ifndef IFMRMCMC_H
#define IFMRMCMC_H

#include <array>
#include <string>
#include <fstream>

#include "evolve.hpp"

const double MIN_MASS1    = 0.15;
const int N_MS_MASS1      = 500;
const int N_MS_MASS_RATIO = 20;
const int N_WD_MASS1      = 8000;  /* evaluate white dwarfs on a finer grid? */
const int CLUS_FILE = 10;
const int ALLOC_CHUNK = 5;
const int nSave = 10;             /*changed from 100 to 10 */

struct ifmrMcmcControl
{
    std::ifstream rData;
    std::ofstream resFile;
    std::ofstream burninFile;
    std::string clusterFilename;

    int nIter;                  // number of post burn-in iterations
    int burnIter;               // total number of burn-in iterations
    int thin;
    std::array<double, NPARAMS> priorVar; // Rewrite the way this is handled... This is ugly.
    int numFilts;
    double propMatrix[NPARAMS][NPARAMS];
};

#endif
