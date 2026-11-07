// Copyright (c) Embry-Riddle Aeronautical University. All rights reserved.

#ifndef EEPPOWERCOMPARE_SETTINGS_H
#define EEPPOWERCOMPARE_SETTINGS_H

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
    int minPower = -3;
    int maxPower = 3;
    int referencePower = 0;
    int gridPoints = 10;
    unsigned int threads = 4;
    bool noBinaries = false;

  private:
    void printUsage();
};

#endif
