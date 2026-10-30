#ifndef SINGLEPOP_SETTINGS_H
#define SINGLEPOP_SETTINGS_H

#include <array>
#include <string>

#include "constants.hpp"
#include "GlobalSettings.hpp"

class Settings
{
  public:
    void loadSettings(int argc, char** argv, const std::string& defaultFile = "base9.yaml");

    void fromYaml(const std::string&);
    void fromCLI(int, char**);

    GlobalSettings global;

    // App-specific flags
    bool noBinaries = false;
    bool onlyWDs = false;
    int eepInterpolationPower = 0;
    int wdInterpolationPower = 0;

    struct SinglePopMcmcSettings
    {
        int burnIter;
        int maxIter;
        int thin;

        int stage3Iter;

        int adaptiveBigSteps;
        int trialIter;

        bool bigStepBurnin = false;

        std::array<double, NPARAMS> stepSize;
    };

    struct SinglePopMcmcSettings singlePopMcmc;

  private:
    void printUsage();
};

#endif
