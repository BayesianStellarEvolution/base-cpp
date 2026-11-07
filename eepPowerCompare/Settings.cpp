// Copyright (c) Embry-Riddle Aeronautical University. All rights reserved.

#include <string>
#include <iostream>

#include <cstring>
#include <getopt.h>

#include "Settings.hpp"

using std::cerr;
using std::endl;
using std::istringstream;
using std::string;

// Application-specific options (file scope for access from loadSettings)
static struct option long_options[] = {
    {"minPower", required_argument, 0, 0xA1},
    {"maxPower", required_argument, 0, 0xA2},
    {"referencePower", required_argument, 0, 0xA3},
    {"gridPoints", required_argument, 0, 0xA4},
    {"threads", required_argument, 0, 0xA5},
    {"noBinaries", no_argument, 0, 0xA6},
    {"help", no_argument, 0, 0xDB},
    {0, 0, 0, 0}
};


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
        GlobalSettings::printUsage(wantsVerbose);
        printUsage();
        exit(EXIT_SUCCESS);
    }

    global.fromCLI(argc, argv, long_options);

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

    global.fromCLI(argc, argv, long_options);
    fromCLI(argc, argv);
}


void Settings::fromYaml(const string& yamlFile)
{
    // No YAML config for this app yet - could add later if needed
}


void Settings::fromCLI(int argc, char **argv)
{
    char **t_argv = new char*[argc];

    for (int i = 0; i < argc; i++)
    {
        t_argv[i] = new char[strlen(argv[i]) + 1];
        strcpy(t_argv[i], argv[i]);
    }

    int c, option_index;

    optind = 0;
    opterr = 0;  // Suppress errors - GlobalSettings handles unknown option reporting

    while ((c = getopt_long(argc, t_argv, "", long_options, &option_index)) != (-1))
    {
        switch (c)
        {
            case 0xA1:
                istringstream(string(optarg)) >> minPower;
                break;

            case 0xA2:
                istringstream(string(optarg)) >> maxPower;
                break;

            case 0xA3:
                istringstream(string(optarg)) >> referencePower;
                break;

            case 0xA4:
                istringstream(string(optarg)) >> gridPoints;
                break;

            case 0xA5:
                istringstream(string(optarg)) >> threads;
                if (!GlobalSettings::validateThreads(threads))
                {
                    cerr << "You must have at least one thread" << endl;
                    exit(EXIT_FAILURE);
                }
                break;

            case 0xA6:
                noBinaries = true;
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
    cerr << "\neepPowerCompare - EEP Interpolation Power Comparison\n";
    cerr << "====================================================\n\n";

    cerr << "Compares log-likelihood values across different EEP interpolation\n";
    cerr << "powers to assess whether lower powers preserve the likelihood surface.\n\n";

    cerr << "Application-Specific Options:\n";
    cerr << "  --minPower <int>        Minimum EEP power to test (default: -3)\n";
    cerr << "  --maxPower <int>        Maximum EEP power to test (default: 3)\n";
    cerr << "  --referencePower <int>  Reference power for comparison (default: 0)\n";
    cerr << "  --gridPoints <int>      Grid points per parameter (default: 5)\n";
    cerr << "  --threads <int>         Number of threads (default: 4)\n";
    cerr << "  --noBinaries            Use no-binary likelihood\n";
    cerr << "\n";
}
