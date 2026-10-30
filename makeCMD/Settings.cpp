#include <string>
#include <iostream>

#include <cstring>
#include <getopt.h>

#include "Settings.hpp"

using std::cerr;
using std::endl;
using std::istringstream;
using std::string;
using YAML::Node;
using YAML::LoadFile;


void Settings::loadSettings(int argc, char** argv, const string& defaultFile)
{
    // Check for --help before global processing (which would exit)
    bool wantsHelp = false;
    bool wantsVerbose = false;
    for (int i = 1; i < argc; ++i)
    {
        if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0)
            wantsHelp = true;
        if (strcmp(argv[i], "--verbose") == 0)
            wantsVerbose = true;
    }
    if (wantsHelp)
    {
        GlobalSettings::printUsage();
        printUsage();
        if (wantsVerbose)
            GlobalSettings::printVerboseUsage();
        exit(EXIT_SUCCESS);
    }

    global.fromCLI(argc, argv);

    if (!global.files.config.empty())
    {
        global.fromYaml(global.files.config);
        fromYaml(global.files.config);
    }
    else
    {
        global.fromYaml(defaultFile);
        fromYaml(defaultFile);
    }

    global.fromCLI(argc, argv);
    fromCLI(argc, argv);
}


void Settings::fromYaml(const string& yamlFile)
{
    // makeCMD doesn't have its own YAML section, but we check for simCluster.noWDs
    // for compatibility with shared config files
    Node configNode;

    try
    {
        configNode = LoadFile(yamlFile);
    }
    catch (YAML::BadFile &e)
    {
        // Config already loaded by global, just return
        return;
    }

    // Optional: load noWDs from simCluster section if present
    if (configNode["simCluster"])
    {
        Node simConfNode = configNode["simCluster"];
        if (simConfNode["noWDs"])
        {
            simCluster.noWDs = simConfNode["noWDs"].as<bool>();
        }
    }
}


void Settings::fromCLI(int argc, char **argv)
{
    char **t_argv = new char*[argc];

    for (int i = 0; i < argc; i++)
    {
        t_argv[i] = new char[strlen(argv[i]) + 1];
        strcpy(t_argv[i], argv[i]);
    }

    static struct option long_options[] = {
        {"noWDs", no_argument, 0, 0xAC},
        {"onlyWDs", no_argument, 0, 0x9F},
        {"eepInterpolationPower", required_argument, 0, 0xBC},
        {"wdInterpolationPower", required_argument, 0, 0xBB},
        {"help", no_argument, 0, 0xDB},
        {0, 0, 0, 0}
    };

    int c, option_index;

    optind = 0;

    while ((c = getopt_long(argc, t_argv, "", long_options, &option_index)) != (-1))
    {
        switch (c)
        {
            case 0xAC:
                simCluster.noWDs = true;
                break;

            case 0x9F:
                onlyWDs = true;
                break;

            case 0xBC:
                istringstream(string(optarg)) >> eepInterpolationPower;
                if (!GlobalSettings::validateEepPower(eepInterpolationPower))
                {
                    cerr << "EEP interpolation power must be " << GlobalSettings::EEP_POWER_MIN << " or greater." << endl;
                    exit(EXIT_FAILURE);
                }
                break;

            case 0xBB:
                istringstream(string(optarg)) >> wdInterpolationPower;
                if (!GlobalSettings::validateWdPower(wdInterpolationPower))
                {
                    cerr << "WD interpolation power must be " << GlobalSettings::WD_POWER_MIN << " or greater." << endl;
                    exit(EXIT_FAILURE);
                }
                break;

            case 0xDB:
                printUsage();
                exit(EXIT_SUCCESS);

            case '?':
                // Ignore unknown options (handled by global settings)
                break;

            default:
                break;
        }
    }

    for (int i = 0; i < argc; i++)
    {
        delete[] t_argv[i];
    }
    delete[] t_argv;
}


void Settings::printUsage()
{
    cerr << "\nmakeCMD / makeCMDs / makeCMDwd - Generate Color-Magnitude Diagrams\n";
    cerr << "====================================================================\n\n";

    cerr << "Generate theoretical color-magnitude diagrams from isochrones.\n";
    cerr << "  makeCMD   - Single CMD output\n";
    cerr << "  makeCMDs  - CMD with sampled cluster parameters\n";
    cerr << "  makeCMDwd - WD-only CMD\n\n";

    cerr << "Application-Specific Options:\n";
    cerr << "  --noWDs                   Do not include white dwarfs in CMD\n";
    cerr << "  --onlyWDs                 Integrate only above AGB tip (WD-only runs)\n";
    cerr << "  --eepInterpolationPower <int>  EEP interpolation resolution\n";
    cerr << "  --wdInterpolationPower <int>   WD mass grid resolution\n";
    cerr << "\n";

}
