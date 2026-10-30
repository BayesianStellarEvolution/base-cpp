#ifndef SAMPLEWDMASS_SETTINGS_H
#define SAMPLEWDMASS_SETTINGS_H

#include <string>

#include "GlobalSettings.hpp"

class Settings
{
  public:
    void loadSettings(int argc, char** argv, const std::string& defaultFile = "base9.yaml");

    void fromCLI(int, char**);

    GlobalSettings global;

    // App-specific flags
    bool onlyWDs = false;

  private:
    void printUsage();
};

#endif
