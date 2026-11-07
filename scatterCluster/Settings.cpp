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
using std::pair;
using YAML::Node;
using YAML::LoadFile;

// Application-specific options (file scope for access from loadSettings)
static struct option long_options[] = {
    {"scatterConfig", required_argument, 0, 'c'},
    {"output", required_argument, 0, 'o'},
    {"brightLimit", required_argument, 0, 0xE4},
    {"faintLimit", required_argument, 0, 0xE3},
    {"relevantFilt", required_argument, 0, 0xE2},
    {"limitS2N", required_argument, 0, 0xE1},
    {"help", no_argument, 0, 'h'},
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

    // Require --scatterConfig
    if (scatterCluster.scatterConfig.empty())
    {
        cerr << "Error: --scatterConfig is required\n" << endl;
        printUsage();
        exit(EXIT_FAILURE);
    }

    // Load scatter configuration
    loadScatterConfig(scatterCluster.scatterConfig);
}


void Settings::fromYaml(const string& yamlFile)
{
    // No scatterCluster-specific settings in base9.yaml anymore
    // All scatter settings come from the separate scatter config file
}


void Settings::loadScatterConfig(const string& configFile)
{
    Node configNode;

    try
    {
        configNode = LoadFile(configFile);
    }
    catch (YAML::BadFile &e)
    {
        cerr << "Scatter configuration file '" << configFile << "' not found." << endl;
        exit(1);
    }

    scatterCluster.brightLimit = global.getOrRequest<double>(configNode, "brightLimit");
    scatterCluster.faintLimit = global.getOrRequest<double>(configNode, "faintLimit");
    scatterCluster.relevantFilt = global.getOrRequest<int>(configNode, "relevantFilt");
    scatterCluster.limitS2N = global.getOrRequest<double>(configNode, "limitS2N");
    scatterCluster.minSigma = global.getOrRequest<double>(configNode, "minSigma");

    {
        auto tNode = global.getNode(configNode, "exposures");
        scatterCluster.exposures = tNode.as<map<string, double>>();
    }

    {
        auto tNode = global.getNode(configNode, "s2nCoeffs");
        for (auto it = tNode.begin(); it != tNode.end(); ++it)
        {
            string filter = it->first.as<string>();
            auto coeffs = it->second.as<std::vector<double>>();
            if (coeffs.size() >= 2)
            {
                scatterCluster.s2nCoeffs[filter] = {coeffs[0], coeffs[1]};
            }
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

    int c, option_index;

    optind = 0;

    while ((c = getopt_long(argc, t_argv, "o:", long_options, &option_index)) != (-1))
    {
        switch (c)
        {
            case 'c':
                scatterCluster.scatterConfig = optarg;
                break;

            case 'o':
                scatterCluster.outputFile = optarg;
                break;

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

            case 'h':
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

    cerr << "Usage: scatterCluster --scatterConfig <file> [options]\n\n";

    cerr << "Required:\n";
    cerr << "  --scatterConfig <path>    Scatter specification YAML file\n";
    cerr << "\n";

    cerr << "Optional:\n";
    cerr << "  -o, --output <path>       Output file path (default: <outputFileBase>.scatter.out)\n";
    cerr << "  --photFile <path>         Input photometry file (default: <outputFileBase>.sim.out)\n";
    cerr << "\n";

    cerr << "Overrides (from scatter config):\n";
    cerr << "  --brightLimit <float>     Bright magnitude limit\n";
    cerr << "  --faintLimit <float>      Faint magnitude limit\n";
    cerr << "  --relevantFilt <int>      Filter index for limits (0 = bluest)\n";
    cerr << "  --limitS2N <float>        Signal-to-noise threshold\n";
    cerr << "\n";

    cerr << "Example scatter configs: conf/scatter_uncrowded.yaml, conf/scatter_crowded.yaml\n";
    cerr << "\n";
}
