#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>

#include "yaml-cpp/yaml.h"

struct MainSequenceSettings
{
    int filterSet;
    int msRgbModel;
};

struct WhiteDwarfSettings
{
    int ifmr;
    int wdModel;
    double carbonicity;
    double M_wd_up;
};

struct BrownDwarfSettings
{
    int bdModel;
};

struct MpiMcmcSettings
{
    int burnIter;
    int maxIter;
    int thin;
};

struct SimClusterSettings
{
    int nStars;
    int nFieldStars;
    int nBrownDwarfs;
    int percentBinary;      // Fraction * 100
    int percentDB;          // Fraction * 100
};

struct ScatterClusterSettings
{
    int relevantFilt;
    double exposures[14];
    double brightLimit;
    double faintLimit;
    double limitS2N;
};

struct ClusterSigmas
{
    double Fe_H;
    double distMod;
    double Av;
    double Y;
};

struct ClusterSettings
{
    double Fe_H;
    double distMod;
    double Av;
    double Y;

    struct ClusterSigmas sigma;

    double logClusAge;

    double minMag;
    double maxMag;
    int index;
};

struct Files
{
    std::string phot;
    std::string output;
    std::string scatter;
    std::string config;
    std::string models;
};

struct Settings
{
    int seed;
    int verbose;

    struct Files files;
    struct MainSequenceSettings mainSequence;
    struct WhiteDwarfSettings whiteDwarf;
    struct BrownDwarfSettings brownDwarf;
    struct MpiMcmcSettings mpiMcmc;
    struct ClusterSettings cluster;
    struct SimClusterSettings simCluster;
    struct ScatterClusterSettings scatterCluster;
};

void makeSettings (const std::string, struct Settings &);
void settingsFromCLI (int, char **, struct Settings &);

template <typename T> T getDefault (YAML::Node &, std::string &&, T);
template <typename T> T getOrDie (YAML::Node &, std::string &&);
YAML::Node getNode (YAML::Node &, std::string &&);
[[noreturn]] void exitWith (std::string &&);
#endif
