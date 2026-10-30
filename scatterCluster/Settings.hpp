#ifndef SCATTERCLUSTER_SETTINGS_H
#define SCATTERCLUSTER_SETTINGS_H

#include <map>
#include <string>

#include "GlobalSettings.hpp"

class Settings
{
  public:
    void loadSettings(int argc, char** argv, const std::string& defaultFile = "base9.yaml");

    void fromYaml(const std::string&);
    void fromCLI(int, char**);

    GlobalSettings global;

    struct ScatterClusterSettings
    {
        int relevantFilt;
        std::map<std::string, double> exposures;
        double brightLimit;
        double faintLimit;
        double limitS2N;
        bool crowded = false;
    };

    struct ScatterClusterSettings scatterCluster;

  private:
    void printUsage();
};

#endif
