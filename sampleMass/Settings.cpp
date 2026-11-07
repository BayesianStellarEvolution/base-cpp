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

// Application-specific options (file scope for access from loadSettings)
static struct option long_options[] = {
    {"deltaMass", required_argument, 0, 0xCB},
    {"deltaMassRatio", required_argument, 0, 0xCA},
    {"noBinaries", no_argument, 0, 0xAE},
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

    Node sampleMassNode = global.getNode(configNode, "sampleMass");

    sampleMass.burnIters = global.getOrRequest<int>(sampleMassNode, "burnIters");
    sampleMass.iters = global.getOrRequest<int>(sampleMassNode, "iters");
    sampleMass.deltaMass = global.getOrRequest<double>(sampleMassNode, "deltaMass");
    sampleMass.deltaMassRatio = global.getOrRequest<double>(sampleMassNode, "deltaMassRatio");
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
            case 0xCB:
                istringstream(string(optarg)) >> sampleMass.deltaMass;
                break;

            case 0xCA:
                istringstream(string(optarg)) >> sampleMass.deltaMassRatio;
                break;

            case 0xAE:
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
    cerr << "\nsampleMass - Sample Stellar Masses\n";
    cerr << "===================================\n\n";

    cerr << "Samples stellar masses from posterior distributions given\n";
    cerr << "cluster parameters and observed photometry.\n\n";

    cerr << "Application-Specific Options:\n";
    cerr << "  --deltaMass <float>       Step size for primary mass sampling\n";
    cerr << "  --deltaMassRatio <float>  Step size for mass ratio sampling\n";
    cerr << "  --noBinaries              Disable binary star integration\n";
    cerr << "\n";

    cerr << "Additional Parameters (from YAML config):\n";
    cerr << "  burnIters                 Number of burn-in iterations\n";
    cerr << "  iters                     Number of sampling iterations\n";
    cerr << "\n";

}
