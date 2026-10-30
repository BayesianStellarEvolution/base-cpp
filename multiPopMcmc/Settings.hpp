#ifndef MULTIPOP_SETTINGS_H
#define MULTIPOP_SETTINGS_H

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
    int eepInterpolationPower = 0;
    int wdInterpolationPower = 0;
    unsigned int threads = std::numeric_limits<int>::max();

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

    struct MultiPopMcmcSettings
    {
        double YA_start;
        double YB_start;
        double lambda_start;

        double YA_lo;
        double YA_hi;
        double YB_hi;

        double lambdaStep;
    };

    struct SinglePopMcmcSettings singlePopMcmc;
    struct MultiPopMcmcSettings multiPopMcmc;

  private:
    void printUsage();
};

#endif
