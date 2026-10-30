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
    Node configNode;

    try
    {
        configNode = LoadFile(yamlFile);
    }
    catch (YAML::BadFile &e)
    {
        cerr << "Configuration file '" << yamlFile << "' not found." << endl;
        exit(1);
    }

    Node simConfNode = global.getNode(configNode, "simCluster");

    simCluster.nStars = global.getOrRequest<int>(simConfNode, "nStars");
    simCluster.percentBinary = global.getOrRequest<int>(simConfNode, "percentBinary");
    simCluster.percentDB = global.getOrRequest<int>(simConfNode, "percentDB");
    simCluster.nFieldStars = global.getOrRequest<int>(simConfNode, "nFieldStars");
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
        {"nStars", required_argument, 0, 0xE9},
        {"percentBinary", required_argument, 0, 0xE8},
        {"percentDB", required_argument, 0, 0xE7},
        {"nFieldStars", required_argument, 0, 0xE6},
        {"noWDs", no_argument, 0, 0xAC},
        {"help", no_argument, 0, 0xDB},
        {0, 0, 0, 0}
    };

    int c, option_index;

    optind = 0;

    while ((c = getopt_long(argc, t_argv, "", long_options, &option_index)) != (-1))
    {
        switch (c)
        {
            case 0xE9:
                istringstream(string(optarg)) >> simCluster.nStars;
                break;

            case 0xE8:
                istringstream(string(optarg)) >> simCluster.percentBinary;
                break;

            case 0xE7:
                istringstream(string(optarg)) >> simCluster.percentDB;
                break;

            case 0xE6:
                istringstream(string(optarg)) >> simCluster.nFieldStars;
                break;

            case 0xAC:
                simCluster.noWDs = true;
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
    cerr << "\nsimCluster - Simulate Stellar Cluster\n";
    cerr << "======================================\n\n";

    cerr << "Generates synthetic photometry for a simulated stellar cluster\n";
    cerr << "based on cluster parameters and stellar evolution models.\n\n";

    cerr << "Application-Specific Options:\n";
    cerr << "  --nStars <int>            Number of stellar systems to simulate\n";
    cerr << "  --percentBinary <int>     Binary fraction (0-100)\n";
    cerr << "  --percentDB <int>         Fraction of WDs with He atmosphere (0-100)\n";
    cerr << "  --nFieldStars <int>       Number of field stars to add\n";
    cerr << "  --noWDs                   Do not generate white dwarfs\n";
    cerr << "\n";

}
