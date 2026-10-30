#ifndef SIMCLUSTER_SETTINGS_H
#define SIMCLUSTER_SETTINGS_H

#include <string>

#include "GlobalSettings.hpp"

class Settings
{
  public:
    void loadSettings(int argc, char** argv, const std::string& defaultFile = "base9.yaml");

    void fromYaml(const std::string&);
    void fromCLI(int, char**);

    GlobalSettings global;

    struct SimClusterSettings
    {
        bool noWDs = false;

        int nStars;
        int nFieldStars;
        int percentBinary;      // Fraction * 100
        int percentDB;          // Fraction * 100
    };

    struct SimClusterSettings simCluster;

  private:
    void printUsage();
};

#endif
