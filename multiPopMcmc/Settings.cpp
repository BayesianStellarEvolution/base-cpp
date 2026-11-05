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

    Node singlePopConfNode = global.getNode(configNode, "singlePopMcmc");
    Node mpiAdaptiveNode = global.getNode(singlePopConfNode, "adaptive");
    Node mpiStepNode = global.getNode(singlePopConfNode, "stepSizes");
    Node multiPopConfNode = global.getNode(configNode, "multiPopMcmc");

    // Load singlePopMcmc settings (shared step sizes)
    singlePopMcmc.burnIter = global.getOrRequest<int>(singlePopConfNode, "stage2IterMax");
    singlePopMcmc.stage3Iter = global.getOrRequest<int>(singlePopConfNode, "stage3Iter");
    singlePopMcmc.maxIter = global.getOrRequest<int>(singlePopConfNode, "runIter");
    singlePopMcmc.thin = global.getOrRequest<int>(singlePopConfNode, "thin");

    singlePopMcmc.adaptiveBigSteps = global.getOrRequest<int>(mpiAdaptiveNode, "bigStepIter");
    singlePopMcmc.trialIter = global.getOrRequest<int>(mpiAdaptiveNode, "trialIter");

    if (singlePopMcmc.trialIter <= 0)
        global.exitWith("singlePopMcmc:adaptive:trialIter must be greater than 0");

    singlePopMcmc.stepSize[AGE] = global.getOrRequest<double>(mpiStepNode, "age");
    singlePopMcmc.stepSize[FEH] = global.getOrRequest<double>(mpiStepNode, "Fe_H");
    singlePopMcmc.stepSize[MOD] = global.getOrRequest<double>(mpiStepNode, "distMod");
    singlePopMcmc.stepSize[ABS] = global.getOrRequest<double>(mpiStepNode, "Av");
    singlePopMcmc.stepSize[YYY] = global.getOrRequest<double>(mpiStepNode, "Y");
    singlePopMcmc.stepSize[CARBONICITY] = global.getOrRequest<double>(mpiStepNode, "carbonicity");
    singlePopMcmc.stepSize[IFMR_INTERCEPT] = global.getOrRequest<double>(mpiStepNode, "ifmrIntercept");
    singlePopMcmc.stepSize[IFMR_SLOPE] = global.getOrRequest<double>(mpiStepNode, "ifmrSlope");
    singlePopMcmc.stepSize[IFMR_QUADCOEF] = global.getOrRequest<double>(mpiStepNode, "ifmrQuadCoef");

    // Load multiPopMcmc-specific settings
    multiPopMcmc.YA_start = global.getOrRequest<double>(multiPopConfNode, "YA_start");
    multiPopMcmc.YB_start = global.getOrRequest<double>(multiPopConfNode, "YB_start");
    multiPopMcmc.lambda_start = global.getOrRequest<double>(multiPopConfNode, "lambda_start");

    multiPopMcmc.YA_lo = global.getOrRequest<double>(multiPopConfNode, "YA_lo");
    multiPopMcmc.YA_hi = global.getOrRequest<double>(multiPopConfNode, "YA_hi");
    multiPopMcmc.YB_hi = global.getOrRequest<double>(multiPopConfNode, "YB_hi");

    multiPopMcmc.lambdaStep = global.getOrRequest<double>(multiPopConfNode, "lambdaStep");
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
        {"burnIter", required_argument, 0, 0xEC},
        {"maxIter", required_argument, 0, 0xEB},
        {"thin", required_argument, 0, 0xEA},
        {"bigStepBurnin", no_argument, 0, 0xCF},
        {"noBinaries", no_argument, 0, 0xAE},
        {"overrideBounds", no_argument, 0, 0xAD},
        {"eepInterpolationPower", required_argument, 0, 0xBC},
        {"wdInterpolationPower", required_argument, 0, 0xBB},
        {"threads", required_argument, 0, 0xCE},
        {"help", no_argument, 0, 0xDB},
        {0, 0, 0, 0}
    };

    int c, option_index;

    optind = 0;

    while ((c = getopt_long(argc, t_argv, "", long_options, &option_index)) != (-1))
    {
        switch (c)
        {
            case 0xEC:
                istringstream(string(optarg)) >> singlePopMcmc.burnIter;
                break;

            case 0xEB:
                istringstream(string(optarg)) >> singlePopMcmc.maxIter;
                break;

            case 0xEA:
                istringstream(string(optarg)) >> singlePopMcmc.thin;
                break;

            case 0xCF:
                singlePopMcmc.bigStepBurnin = true;
                break;

            case 0xAE:
                noBinaries = true;
                break;

            case 0xAD:
                overrideBounds = true;
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

            case 0xCE:
                istringstream(string(optarg)) >> threads;
                if (!GlobalSettings::validateThreads(threads))
                {
                    cerr << "You must have at least one thread" << endl;
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
    cerr << "\nmultiPopMcmc - Multiple Population MCMC Fitting\n";
    cerr << "================================================\n\n";

    cerr << "Performs Bayesian analysis of stellar cluster photometry using\n";
    cerr << "Markov Chain Monte Carlo methods with a dual stellar population model.\n";
    cerr << "Fits two populations with different helium abundances.\n\n";

    cerr << "Application-Specific Options:\n";
    cerr << "  --burnIter <int>          Maximum Stage 2 burn-in iterations\n";
    cerr << "  --maxIter <int>           Total MCMC iterations to run\n";
    cerr << "  --thin <int>              Thinning factor for output\n";
    cerr << "  --bigStepBurnin           Use large steps during burn-in\n";
    cerr << "  --noBinaries              Disable binary star integration\n";
    cerr << "  --overrideBounds          Skip model boundary validation\n";
    cerr << "  --eepInterpolationPower <int>  EEP interpolation resolution\n";
    cerr << "  --wdInterpolationPower <int>   WD mass grid resolution\n";
    cerr << "  --threads <int>           Number of threads to use\n";
    cerr << "\n";

    cerr << "Multi-Population Parameters (from YAML config):\n";
    cerr << "  YA_start, YB_start        Starting helium for populations A and B\n";
    cerr << "  lambda_start              Starting population fraction\n";
    cerr << "  YA_lo, YA_hi, YB_hi       Helium bounds\n";
    cerr << "  lambdaStep                Lambda step size\n";
    cerr << "\n";

}
