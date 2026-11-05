#ifndef GLOBAL_SETTINGS_H
#define GLOBAL_SETTINGS_H

#include <array>
#include <string>
#include <vector>

#include "constants.hpp"
#include "yaml-cpp/yaml.h"

enum class Backend { File, Sqlite };

class GlobalSettings
{
  public:
    void loadSettings(int argc, char** argv, const std::string& defaultFile = "base9.yaml");

    void fromYaml (const std::string&);
    void fromCLI (int, char **);

    bool modIsParallax = false;
    bool allowInvalidModels = false;

    bool details = false;
    bool verbose = false;
    bool veryVerbose = false;

    bool includeBurnin = false;

    bool ignoreLowSigma = false;     // Disables low sigma warning
    bool allowNegativeSigma = false; // Enables disabling individual photometric observations in the photometry file

    std::vector<std::string> include;
    std::vector<std::string> exclude;

    uint32_t seed = std::numeric_limits<uint32_t>::max();


    struct MainSequenceSettings
    {
        // Path to MS model, relative to "models"
        std::string modelFile = "";
    };

    struct WhiteDwarfSettings
    {
        int ifmr;
        WdModel wdModel;
        WdAtmosphereModelSet wdAtmosphereModel;
        double M_wd_up;
    };

    struct ClusterValues
    {
        double Fe_H;
        double distMod;
        double Av;
        double Y;
        double carbonicity;
        double logAge;
    };

    struct ClusterSettings
    {
        struct ClusterValues priorMeans;
        struct ClusterValues priorSigma;
        struct ClusterValues starting;

        double minMag;
        double maxMag;
        int index;
    };

    struct Files
    {
        Backend backend;

        std::string phot;
        std::string output;
        std::string scatter;
        std::string config;
        std::string models;
    };

    struct Files files;
    struct MainSequenceSettings mainSequence;
    struct WhiteDwarfSettings whiteDwarf;
    struct ClusterSettings cluster;

    // Helper functions for YAML parsing (public for use by app-specific Settings)
    template <typename T> T getDefault (YAML::Node &, std::string, T);
    template <typename T> T getOrDie (YAML::Node &, std::string);
    template <typename T> T getOrRequest (YAML::Node &, std::string);
    YAML::Node getNode (YAML::Node &, std::string);
    [[noreturn]] void exitWith (std::string);

    static void printUsage();
    static void printVerboseUsage();

    // Static helpers for interpolation power validation and step calculation
    static constexpr int EEP_POWER_MIN = -5;
    static constexpr int WD_POWER_MIN = -7;

    static bool validateEepPower(int power) { return power >= EEP_POWER_MIN; }
    static bool validateWdPower(int power) { return power >= WD_POWER_MIN; }
    static bool validateThreads(unsigned int threads) { return threads >= 1; }

    static int calcEepSteps(int eepInterpolationPower) {
        return eepInterpolationPower == 0 ? 80
             : eepInterpolationPower > 0 ? 64 << eepInterpolationPower
                                         : 64 >> -eepInterpolationPower;
    }

    static int calcWdMassSteps(int wdInterpolationPower) {
        return wdInterpolationPower == 0 ? 8000
             : wdInterpolationPower > 0 ? 4096 << wdInterpolationPower
                                        : 4096 >> -wdInterpolationPower;
    }
};

#endif
