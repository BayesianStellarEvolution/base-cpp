#ifndef SCATTERCLUSTER_SETTINGS_H
#define SCATTERCLUSTER_SETTINGS_H

#include <map>
#include <string>
#include <utility>

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
        std::string scatterConfig;  // Required scatter specification YAML
        std::string outputFile;     // Optional output file override

        // Loaded from scatter config
        int relevantFilt = 0;
        double brightLimit = -10;
        double faintLimit = 30;
        double limitS2N = 1;
        double minSigma = 0.01;
        std::map<std::string, double> exposures;
        std::map<std::string, std::pair<double, double>> s2nCoeffs;
    };

    struct ScatterClusterSettings scatterCluster;

    void loadScatterConfig(const std::string& configFile);

  private:
    void printUsage();
};

#endif
