#ifndef SAMPLEMASS_SETTINGS_H
#define SAMPLEMASS_SETTINGS_H

#include <string>

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

    struct SampleMassSettings
    {
        int iters;
        int burnIters;

        double deltaMass;
        double deltaMassRatio;
    };

    struct SampleMassSettings sampleMass;

  private:
    void printUsage();
};

#endif
