#include <string>
#include <iostream>
#include <map>

#include <cstring>
#include <getopt.h>

#include "Settings.hpp"

using std::cerr;
using std::endl;
using std::istringstream;
using std::string;
using std::map;
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

    Node scatterConfNode = global.getNode(configNode, "scatterCluster");

    scatterCluster.brightLimit = global.getOrRequest<double>(scatterConfNode, "brightLimit");
    scatterCluster.faintLimit = global.getOrRequest<double>(scatterConfNode, "faintLimit");
    scatterCluster.relevantFilt = global.getOrRequest<int>(scatterConfNode, "relevantFilt");
    scatterCluster.limitS2N = global.getOrRequest<double>(scatterConfNode, "limitS2N");
    scatterCluster.crowded = global.getOrRequest<bool>(scatterConfNode, "crowded");

    {
        auto tNode = global.getNode(scatterConfNode, "exposures");
        scatterCluster.exposures = tNode.as<map<string, double>>();
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
        {"brightLimit", required_argument, 0, 0xE4},
        {"faintLimit", required_argument, 0, 0xE3},
        {"relevantFilt", required_argument, 0, 0xE2},
        {"limitS2N", required_argument, 0, 0xE1},
        {"scatterFile", required_argument, 0, 0xDE},
        {"help", no_argument, 0, 0xDB},
        {0, 0, 0, 0}
    };

    int c, option_index;

    optind = 0;

    while ((c = getopt_long(argc, t_argv, "", long_options, &option_index)) != (-1))
    {
        switch (c)
        {
            case 0xE4:
                istringstream(string(optarg)) >> scatterCluster.brightLimit;
                break;

            case 0xE3:
                istringstream(string(optarg)) >> scatterCluster.faintLimit;
                break;

            case 0xE2:
                istringstream(string(optarg)) >> scatterCluster.relevantFilt;
                break;

            case 0xE1:
                istringstream(string(optarg)) >> scatterCluster.limitS2N;
                break;

            case 0xDE:
                global.files.scatter = optarg;
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
    cerr << "\nscatterCluster - Apply Photometric Scatter\n";
    cerr << "===========================================\n\n";

    cerr << "Adds realistic photometric scatter to simulated cluster data\n";
    cerr << "based on exposure times and signal-to-noise characteristics.\n\n";

    cerr << "Application-Specific Options:\n";
    cerr << "  --brightLimit <float>     Bright magnitude limit (remove bright stars)\n";
    cerr << "  --faintLimit <float>      Faint magnitude limit (remove faint stars)\n";
    cerr << "  --relevantFilt <int>      Filter index for limits (0 = bluest)\n";
    cerr << "  --limitS2N <float>        Signal-to-noise threshold\n";
    cerr << "  --scatterFile <path>      Output scatter file path\n";
    cerr << "\n";

    cerr << "Additional Parameters (from YAML config):\n";
    cerr << "  exposures                 Exposure times per filter\n";
    cerr << "  crowded                   Use crowded-field S/N model\n";
    cerr << "\n";

}
