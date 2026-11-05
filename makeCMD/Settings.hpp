#ifndef MAKECMD_SETTINGS_H
#define MAKECMD_SETTINGS_H

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
    bool onlyWDs = false;
    bool interpolateCMD = false;
    int eepInterpolationPower = 0;
    int wdInterpolationPower = 0;

    // Re-use SimClusterSettings for noWDs flag
    struct SimClusterSettings
    {
        bool noWDs = false;
    };

    struct SimClusterSettings simCluster;

  private:
    void printUsage();
};

#endif
