#include <string>
#include <iostream>
#include <map>
#include <vector>

#include <cstdio>
#include <cstring>
#include <getopt.h>

#include "Base9Config.h"
#include "yaml-cpp/yaml.h"
#include "Settings.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::istringstream;
using std::string;
using std::vector;
using std::map;
using YAML::Node;
using YAML::LoadFile;

// Forward declaration
static void printUsage ();
static void printVersion ();


void Settings::loadSettings(int argc, char** argv, const string& defaultFile)
{
    fromCLI (argc, argv);

    if (!files.config.empty())
    {
        fromYaml (files.config);
    }
    else
    {
        fromYaml (defaultFile);
    }

    fromCLI (argc, argv);
}


void Settings::fromYaml (const string& yamlFile)
{
    Node configNode;

    try
    {
        configNode = LoadFile (yamlFile);
    }
    catch (YAML::BadFile &e)
    {
        cerr << "Configuration file '" << yamlFile << "' not found." << endl;
        exit(1);
    }

    Node generalNode = getNode (configNode, "general");
    Node filesNode = getNode (generalNode, "files");
    Node mainSequenceNode = getNode (generalNode, "main_sequence");
    Node whiteDwarfNode = getNode (generalNode, "white_dwarfs");
    Node clusterNode = getNode (generalNode, "cluster");
    Node startingNode = getNode(clusterNode, "starting");
    Node priorNode = getNode(clusterNode, "priors");
    Node meansNode = getNode (priorNode, "means");
    Node sigmasNode = getNode (priorNode, "sigmas");
    Node singlePopConfNode = getNode (configNode, "singlePopMcmc");
    Node mpiAdaptiveNode = getNode(singlePopConfNode, "adaptive");
    Node mpiStepNode = getNode (singlePopConfNode, "stepSizes");
    Node multiPopConfNode = getNode (configNode, "multiPopMcmc");
    Node cmdConfNode = getNode (configNode, "makeCMD");
    Node simConfNode = getNode (configNode, "simCluster");
    Node isoConfNode = getNode (configNode, "makeIsochrone");
    Node scatterConfNode = getNode (configNode, "scatterCluster");
    Node sampleMassNode = getNode (configNode, "sampleMass");

    mainSequence.modelFile = getOrRequest <string>(mainSequenceNode, "modelFile");

    whiteDwarf.ifmr = getOrRequest <int>(whiteDwarfNode, "ifmr");
    whiteDwarf.wdModel = static_cast<WdModel>(getOrRequest <int>(whiteDwarfNode, "wdModel"));
    whiteDwarf.wdAtmosphereModel = static_cast<WdAtmosphereModelSet>(getOrRequest <int>(whiteDwarfNode, "wdAtmosphereModel"));
    whiteDwarf.M_wd_up = getOrRequest <double>(whiteDwarfNode, "M_wd_up");

    cluster.starting.Fe_H = getOrRequest <double>(startingNode, "Fe_H");
    cluster.priorMeans.Fe_H = getOrRequest <double>(meansNode, "Fe_H");
    cluster.priorSigma.Fe_H = getOrRequest <double>(sigmasNode, "Fe_H");

    if ( ( startingNode["distMod"] || meansNode["distMod"] || sigmasNode["distMod"] ) &&
         ( startingNode["parallax"] || meansNode["parallax"] || sigmasNode["parallax"] ) )
    {
        exitWith("You may not provide starting, prior mean, or prior sigma values for both parallax and distMod simultaneously.");
    }
    else if (startingNode["distMod"] && meansNode["distMod"] && sigmasNode["distMod"])
    {
        modIsParallax = false;

        cluster.starting.distMod = getOrRequest <double>(startingNode, "distMod");
        cluster.priorMeans.distMod = getOrRequest <double>(meansNode, "distMod");
        cluster.priorSigma.distMod = getOrRequest <double>(sigmasNode, "distMod");
    }
    else if (startingNode["parallax"] && meansNode["parallax"] && sigmasNode["parallax"])
    {
        modIsParallax = true;

        cluster.starting.distMod = getOrRequest <double>(startingNode, "parallax");
        cluster.priorMeans.distMod = getOrRequest <double>(meansNode, "parallax");
        cluster.priorSigma.distMod = getOrRequest <double>(sigmasNode, "parallax");
    }
    else
    {
        exitWith("You must provide starting, prior mean, and prior sigma value for either distMod or parallax");
    }

    cluster.starting.Av = getOrRequest <double>(startingNode, "Av");
    cluster.priorMeans.Av = getOrRequest <double>(meansNode, "Av");
    cluster.priorSigma.Av = getOrRequest <double>(sigmasNode, "Av");

    cluster.starting.Y = getOrRequest <double>(startingNode, "Y");
    cluster.priorMeans.Y = getOrRequest <double>(meansNode, "Y");
    cluster.priorSigma.Y = getOrRequest <double>(sigmasNode, "Y");

    cluster.starting.carbonicity = getOrRequest <double>(startingNode, "carbonicity");
    cluster.priorMeans.carbonicity = getOrRequest <double>(meansNode, "carbonicity");
    cluster.priorSigma.carbonicity = getOrRequest <double>(sigmasNode, "carbonicity");

    cluster.starting.logAge = getOrRequest <double>(startingNode, "logAge");
    cluster.priorMeans.logAge = getOrRequest <double>(meansNode, "logAge");
    cluster.priorSigma.logAge = getOrRequest <double>(sigmasNode, "logAge");

    cluster.minMag = getOrRequest <double>(clusterNode, "minMag");
    cluster.maxMag = getOrRequest <double>(clusterNode, "maxMag");
    cluster.index = getOrRequest <int>(clusterNode, "index");

    singlePopMcmc.burnIter = getOrRequest <int>(singlePopConfNode, "stage2IterMax");
    singlePopMcmc.stage3Iter = getOrRequest <int>(singlePopConfNode, "stage3Iter");
    singlePopMcmc.maxIter = getOrRequest <int>(singlePopConfNode, "runIter");
    singlePopMcmc.thin = getOrRequest <int>(singlePopConfNode, "thin");

    singlePopMcmc.adaptiveBigSteps = getOrRequest <int>(mpiAdaptiveNode, "bigStepIter");
    singlePopMcmc.trialIter = getOrRequest <int>(mpiAdaptiveNode, "trialIter");

    if (singlePopMcmc.trialIter <= 0)
        exitWith("singlePopMcmc:adaptive:trialIter must be greater than 0");

    singlePopMcmc.stepSize[AGE] = getOrRequest <double>(mpiStepNode, "age");
    singlePopMcmc.stepSize[FEH] = getOrRequest <double>(mpiStepNode, "Fe_H");
    singlePopMcmc.stepSize[MOD] = getOrRequest <double>(mpiStepNode, "distMod");
    singlePopMcmc.stepSize[ABS] = getOrRequest <double>(mpiStepNode, "Av");
    singlePopMcmc.stepSize[YYY] = getOrRequest <double>(mpiStepNode, "Y");
    singlePopMcmc.stepSize[CARBONICITY] = getOrRequest <double>(mpiStepNode, "carbonicity");
    singlePopMcmc.stepSize[IFMR_INTERCEPT] = getOrRequest <double>(mpiStepNode, "ifmrIntercept");
    singlePopMcmc.stepSize[IFMR_SLOPE] = getOrRequest <double>(mpiStepNode, "ifmrSlope");
    singlePopMcmc.stepSize[IFMR_QUADCOEF] = getOrRequest <double>(mpiStepNode, "ifmrQuadCoef");

    multiPopMcmc.YA_start = getOrRequest <double>(multiPopConfNode, "YA_start");
    multiPopMcmc.YB_start = getOrRequest <double>(multiPopConfNode, "YB_start");
    multiPopMcmc.lambda_start = getOrRequest <double>(multiPopConfNode, "lambda_start");

    multiPopMcmc.YA_lo = getOrRequest <double>(multiPopConfNode, "YA_lo");
    multiPopMcmc.YA_hi = getOrRequest <double>(multiPopConfNode, "YA_hi");
    multiPopMcmc.YB_hi = getOrRequest <double>(multiPopConfNode, "YB_hi");

    multiPopMcmc.lambdaStep = getOrRequest <double>(multiPopConfNode, "lambdaStep");

    simCluster.nStars = getOrRequest <int>(simConfNode, "nStars");
    simCluster.percentBinary = getOrRequest <int>(simConfNode, "percentBinary");
    simCluster.percentDB = getOrRequest <int>(simConfNode, "percentDB");
    simCluster.nFieldStars = getOrRequest <int>(simConfNode, "nFieldStars");
//    simCluster.nBrownDwarfs = getOrRequest <int>(simConfNode, "nBrownDwarfs");

    scatterCluster.brightLimit = getOrRequest <double>(scatterConfNode, "brightLimit");
    scatterCluster.faintLimit = getOrRequest <double>(scatterConfNode, "faintLimit");
    scatterCluster.relevantFilt = getOrRequest <int>(scatterConfNode, "relevantFilt");
    scatterCluster.limitS2N = getOrRequest <double>(scatterConfNode, "limitS2N");
    scatterCluster.crowded  = getOrRequest <bool>(scatterConfNode, "crowded");

    sampleMass.burnIters = getOrRequest <int>(sampleMassNode, "burnIters");
    sampleMass.iters     = getOrRequest <int>(sampleMassNode, "iters");
    sampleMass.deltaMass = getOrRequest <double>(sampleMassNode, "deltaMass");
    sampleMass.deltaMassRatio = getOrRequest <double>(sampleMassNode, "deltaMassRatio");

    {
        auto tNode = getNode(scatterConfNode, "exposures");
        scatterCluster.exposures = tNode.as<map<string, double>>();
    }

    verbose = getOrRequest <int>(generalNode, "verbose");

    // When we switch to C++11, we can change these to std::string and remove most of the cruft
    files.backend = static_cast<Backend>(getOrRequest <int>(filesNode, "backend"));
    files.phot = getOrRequest <string>(filesNode, "photFile");
    files.output = getOrRequest <string>(filesNode, "outputFileBase");
    files.scatter = getOrRequest <string>(filesNode, "scatterFile");
    files.models = getOrRequest <string>(filesNode, "modelDirectory");
}

void commaSeparated(string toParse, std::vector<std::string> &vec)
{
    std::stringstream ss(toParse);
    std::string lid;

    while (std::getline(ss, lid, ',')) {
        vec.push_back(lid);
    }
}

void Settings::fromCLI (int argc, char **argv)
{
    char **t_argv = new char*[argc];

    for (int i = 0; i<argc; i++)
    {
        t_argv[i] = new char[strlen (argv[i]) + 1];

        strcpy (t_argv[i], argv[i]);
    }

    static struct option long_options[] = {
        // These all have to be parsed
        {"msModelFile", required_argument, 0, 0xFE},
        {"ifmr", required_argument, 0, 0xFD},
        {"wdModel", required_argument, 0, 0xFC},
        {"M_wd_up", required_argument, 0, 0xFA},
        {"bdModel", required_argument, 0, 0xF9},
        {"priorFe_H", required_argument, 0, 0xF8},
        {"sigmaFe_H", required_argument, 0, 0xF7},
        {"priorDistMod", required_argument, 0, 0xF6},
        {"sigmaDistMod", required_argument, 0, 0xF5},
        {"priorAv", required_argument, 0, 0xF4},
        {"sigmaAv", required_argument, 0, 0xF3},
        {"priorY", required_argument, 0, 0xF2},
        {"sigmaY", required_argument, 0, 0xF1},
        {"priorLogAge", required_argument, 0, 0xF0},
        {"minMag", required_argument, 0, 0xEF},
        {"maxMag", required_argument, 0, 0xEE},
        {"index", required_argument, 0, 0xED},
        {"burnIter", required_argument, 0, 0xEC},
        {"maxIter", required_argument, 0, 0xEB},
        {"thin", required_argument, 0, 0xEA},
        {"nStars", required_argument, 0, 0xE9},
        {"percentBinary", required_argument, 0, 0xE8},
        {"percentDB", required_argument, 0, 0xE7},
        {"nFieldStars", required_argument, 0, 0xE6},
        {"modelDirectory", required_argument, 0, 0xE5},
        {"brightLimit", required_argument, 0, 0xE4},
        {"faintLimit", required_argument, 0, 0xE3},
        {"relevantFilt", required_argument, 0, 0xE2},
        {"limitS2N", required_argument, 0, 0xE1},
        {"seed", required_argument, 0, 0xE0},
        {"photFile", required_argument, 0, 0xDF},
        {"scatterFile", required_argument, 0, 0xDE},
        {"outputFileBase", required_argument, 0, 0xDD},
        {"config", required_argument, 0, 0xDC},
        {"help", no_argument, 0, 0xDB},
        {"version", no_argument, 0, 0xDA},
        {"bigStepBurnin", no_argument, 0, 0xCF},
        {"threads", required_argument, 0, 0xCE},
        {"priorCarbonicity", required_argument, 0, 0xCD},
        {"sigmaCarbonicity", required_argument, 0, 0xCC},
        {"deltaMass", required_argument, 0, 0xCB},
        {"deltaMassRatio", required_argument, 0, 0xCA},
        {"sigmaLogAge", required_argument, 0, 0xC9},
        {"startingFe_H", required_argument, 0, 0xC8},
        {"startingDistMod", required_argument, 0, 0xC7},
        {"startingAv", required_argument, 0, 0xC6},
        {"startingY", required_argument, 0, 0xC5},
        {"startingLogAge", required_argument, 0, 0xC4},
        {"startingCarbonicity", required_argument, 0, 0xC3},
        {"backend", required_argument, 0, 0xC2},

        {"startingParallax", required_argument, 0, 0xC0},
        {"priorParallax", required_argument, 0, 0xBF},
        {"sigmaParallax", required_argument, 0, 0xBE},

        {"wdAtmosphereModel", required_argument, 0, 0xBD},

        // Various flags
        // These are now handled the same way as the parameters due to occasional compiler weirdness
        {"verbose", no_argument, 0, 0xAF},
        {"noBinaries", no_argument, 0, 0xAE},
        {"overrideBounds", no_argument, 0, 0xAD},
        {"noWDs", no_argument, 0, 0xAC},
        {"exitAfterLoad", no_argument, 0, 0xAB},
        {"details", no_argument, 0, 0xAA},
        {"onlyWDs", no_argument, 0, 0x9F},
        {"allowInvalidModels", no_argument, 0, 0x9E},

        {"exitAfterLoad", no_argument, 0, 0x9D},
        {"ignoreLowSigma", no_argument, 0, 0x9C},
        {"allowNegativeSigma", no_argument, 0, 0x9B},

        {"eepInterpolationPower", required_argument, 0, 0xBC},
        {"wdInterpolationPower", required_argument, 0, 0xBB},
        {"includeBurnin", no_argument, 0, 0xBA},
        {"veryVerbose", no_argument, 0, 0xB9},
        {"stopAfterBurnin", no_argument, 0, 0xB8},
        {"startWithBurnin", required_argument, 0, 0xB7},

        {"include", required_argument, 0, 0x9A},
        {"exclude", required_argument, 0, 0x99},
        {0, 0, 0, 0}
    };

    char modBits = 0;

    int c, option_index;

    optind = 0;

    while ((c = getopt_long (argc, t_argv, "", long_options, &option_index)) != (-1))
    {
        int i;

        switch (c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;

                cout << "option " << long_options[option_index].name;

                if (optarg)
                    cout << " with arg " << optarg;

                cout << endl;
                break;

            case 0xFE:
                mainSequence.modelFile = optarg;
                break;

            case 0xFD:
                istringstream (string (optarg)) >> whiteDwarf.ifmr;
                break;

            case 0xFC:
                istringstream (string (optarg)) >> i;
                whiteDwarf.wdModel = static_cast<WdModel>(i);
                break;

            case 0xBD:
                istringstream (string (optarg)) >> i;
                whiteDwarf.wdAtmosphereModel = static_cast<WdAtmosphereModelSet>(i);
                break;

            case 0xC2:
                istringstream (string (optarg)) >> i;
                files.backend = static_cast<Backend>(i);
                break;

            case 0xFA:
                istringstream (string (optarg)) >> whiteDwarf.M_wd_up;
                break;

            case 0xF8:
                istringstream (string (optarg)) >> cluster.priorMeans.Fe_H;
                break;

            case 0xF7:
                istringstream (string (optarg)) >> cluster.priorSigma.Fe_H;
                break;

            case 0xF4:
                istringstream (string (optarg)) >> cluster.priorMeans.Av;
                break;

            case 0xF3:
                istringstream (string (optarg)) >> cluster.priorSigma.Av;
                break;

            case 0xF2:
                istringstream (string (optarg)) >> cluster.priorMeans.Y;
                break;

            case 0xF1:
                istringstream (string (optarg)) >> cluster.priorSigma.Y;
                break;

            case 0xCD:
                istringstream (string (optarg)) >> cluster.priorMeans.carbonicity;
                break;

            case 0xCC:
                istringstream (string (optarg)) >> cluster.priorSigma.carbonicity;
                break;

            case 0xF0:
                istringstream (string (optarg)) >> cluster.priorMeans.logAge;
                break;

            case 0xEF:
                istringstream (string (optarg)) >> cluster.minMag;
                break;

            case 0xEE:
                istringstream (string (optarg)) >> cluster.maxMag;
                break;

            case 0xED:
                istringstream (string (optarg)) >> cluster.index;
                break;

            case 0xEC:
                istringstream (string (optarg)) >> singlePopMcmc.burnIter;
                break;

            case 0xEB:
                istringstream (string (optarg)) >> singlePopMcmc.maxIter;
                break;

            case 0xEA:
                istringstream (string (optarg)) >> singlePopMcmc.thin;
                break;

            case 0xE9:
                istringstream (string (optarg)) >> simCluster.nStars;
                break;

            case 0xE8:
                istringstream (string (optarg)) >> simCluster.percentBinary;
                break;

            case 0xE7:
                istringstream (string (optarg)) >> simCluster.percentDB;
                break;

            case 0xE6:
                istringstream (string (optarg)) >> simCluster.nFieldStars;
                break;

            case 0xE5:
                files.models = optarg;
                break;

            case 0xE4:
                istringstream (string (optarg)) >> scatterCluster.brightLimit;
                break;

            case 0xE3:
                istringstream (string (optarg)) >> scatterCluster.faintLimit;
                break;

            case 0xE2:
                istringstream (string (optarg)) >> scatterCluster.relevantFilt;
                break;

            case 0xE1:
                istringstream (string (optarg)) >> scatterCluster.limitS2N;
                break;

            case 0xE0:
                istringstream (string (optarg)) >> seed;
                break;

            case 0xDF:
                files.phot = optarg;
                break;

            case 0xDE:
                files.scatter = optarg;
                break;

            case 0xDD:
                files.output = optarg;
                break;

            case 0xDC:
                files.config = optarg;
                break;

            case 0xDB:                  // --help
                printUsage ();
                exit (EXIT_SUCCESS);

            case 0xDA:
                printVersion ();
                exit (EXIT_SUCCESS);

            case 0xCF:
                singlePopMcmc.bigStepBurnin = true;
                break;

            case 0xCE:
                istringstream (string (optarg)) >> threads;
                if (threads <= 0)
                {
                    cerr << "You must have at least one thread" << endl;
                    exit(EXIT_FAILURE);
                }
                break;

            case 0xCB:
                istringstream (string (optarg)) >> sampleMass.deltaMass;
                break;

            case 0xCA:
                istringstream (string (optarg)) >> sampleMass.deltaMassRatio;
                break;

            case 0xC9:
                istringstream (string (optarg)) >> cluster.priorSigma.logAge;
                break;

            case 0xC8:
                istringstream (string (optarg)) >> cluster.starting.Fe_H;
                break;

            case 0xC6:
                istringstream (string (optarg)) >> cluster.starting.Av;
                break;

            case 0xC5:
                istringstream (string (optarg)) >> cluster.starting.Y;
                break;

            case 0xC4:
                istringstream (string (optarg)) >> cluster.starting.logAge;
                break;

            case 0xC3:
                istringstream (string (optarg)) >> cluster.starting.carbonicity;
                break;

            case 0xAF:
                verbose = true;
                break;

            case 0xB9:
                veryVerbose = true;
                break;

            case 0xAE:
                noBinaries = true;
                break;

            case 0xAD:
                overrideBounds = true;
                break;

            case 0xAC:
                simCluster.noWDs = true;
                break;

            case 0xAB:
                exitAfterLoad = true;
                break;

            case 0xAA:
                details = true;
                break;

            case 0x9F:
                onlyWDs = true;
                break;

            case 0x9E:
                allowInvalidModels = true;
                break;

            case 0x9D:
                exitAfterLoad = true;
                break;

            case 0x9C:
                ignoreLowSigma = true;
                break;

            case 0x9B:
                allowNegativeSigma = true;
                break;

            // distMod and parallax have to be handled with special care.
            case 0xC7:
                istringstream (string (optarg)) >> cluster.starting.distMod;
                modBits |= 1;
                break;

            case 0xF6:
                istringstream (string (optarg)) >> cluster.priorMeans.distMod;
                modBits |= 2;
                break;

            case 0xF5:
                istringstream (string (optarg)) >> cluster.priorSigma.distMod;
                modBits |= 4;
                break;

            case 0xC0:
                istringstream (string (optarg)) >> cluster.starting.distMod;
                modBits |= 8;
                break;

            case 0xBF:
                istringstream (string (optarg)) >> cluster.priorMeans.distMod;
                modBits |= 16;
                break;

            case 0xBE:
                istringstream (string (optarg)) >> cluster.priorSigma.distMod;
                modBits |= 32;
                break;

            case 0xBC:
                istringstream (string (optarg)) >> eepInterpolationPower;

                if (eepInterpolationPower < -6)
                {
                    cerr << "EEP interpolation power must be -6 or greater." << endl;
                    exit(EXIT_FAILURE);
                }

                break;

            case 0xBB:
                istringstream (string (optarg)) >> wdInterpolationPower;

                if (wdInterpolationPower < -13)
                {
                    cerr << "WD interpolation power must be -13 or greater." << endl;
                    exit(EXIT_FAILURE);
                }

                break;

            case 0xBA:
                includeBurnin = true;
                break;

            case 0xB8:
                stopAfterBurnin = true;
                break;

            case 0xB7:
                istringstream (string (optarg)) >> startWithBurnin;
                break;

            case 0x9A:
                commaSeparated(string(optarg), include);

                break;

            case 0x99:
                commaSeparated(string(optarg), exclude);

                break;

            case '?':
                // getopt_long already printed an error message.
                printUsage ();
                exit (EXIT_FAILURE);

            default:
                abort ();
        }
    }

    for (int i = 0; i<argc; i++)
    {
        delete[] t_argv[i];
    }
    delete[] t_argv;

    // Print any remaining command line arguments (not options). This is mainly for debugging purposes.
    if (optind < argc)
    {
        cerr << "Unrecognized options: ";
        while (optind<argc)
            cerr << argv[optind++];
        cerr << endl;
    }

    if (modBits == 0)  { return; }
    if (modBits == 7)  { modIsParallax = false; return; }
    if (modBits == 56) { modIsParallax = true;  return; }

    exitWith("You must specify starting, prior mean, and prior sigma values for either parallax or distMod.");
}

template<typename T>
T Settings::getDefault (Node & n, string f, T def)
{
    if (n[f])
    {
        return n[f].as<T> ();
    }
    else
    {
        return def;
    }
}

template<typename T>
T Settings::getOrDie (Node & n, string f)
{
    if (n[f])
    {
        return n[f].as<T> ();
    }
    else
    {
        exitWith ("Field '" + f + "' was not set");
    }
}

template<typename T>
T Settings::getOrRequest (Node & n, string f)
{
    // In the event that there is an issue with loading the YAML file, we
    // request more information from the user.
    //
    // There is no sensible way to request information from the user in some
    // circumstances (e.g., supercomputing clusters, where the client may be run
    // many nodes at once and where the user may not have access to stdin). For
    // these cases, we provide a build flag that makes this function call
    // getOrDie instead.
    #ifndef BUILD_NONINTERACTIVE

    if (n[f])
    {
        return n[f].as<T> ();
    }
    else
    {
        cerr << "No value was found in the YAML file for parameter '" + f + "'.\n Please enter your desired value for this parameter: ";

        string input = "";
        getline(std::cin, input);

        T t;
        istringstream(input) >> t;

        return t;
    }

    #else

    return getOrDie<T>(n, f);

    #endif
}

Node Settings::getNode (Node & n, string f)
{
    if (n[f])
    {
        return n[f];
    }
    else
    {
        exitWith ("Node '" + f + "' was not present\nIs your YAML file up to date?\n");
    }
}

[[noreturn]] void Settings::exitWith (string s)
{
    cerr << s << endl;
    abort ();
}

static void printUsage ()
{
    cerr << "\nUsage:" << endl;
    cerr <<   "======" << endl;

    cerr << "\t--help" << endl;
    cerr << "\t\tPrints help\n" << endl;

    cerr << "\t--version" << endl;
    cerr << "\t\tPrints version string\n" << endl;

    cerr << "\t--config" << endl;
    cerr << "\t\tYAML configuration file\n" << endl;

    cerr << "\t--seed" << endl;
    cerr << "\t\tinitialize the random number generator\n" << endl;

    cerr << "\t--msModelFile" << endl;
    cerr << "\t\tSpecify a model file relative to the --modelDirectory" << endl;

    cerr << "\t--ifmr" << endl;
    cerr << "\t\t0 = Weidemann" << endl;
    cerr << "\t\t1 = Williams" << endl;
    cerr << "\t\t2 = Salaris lin" << endl;
    cerr << "\t\t3 = Salaris pw lin" << endl;
    cerr << "\t\t4 - 11 = tunable" << endl;
    cerr << "\t\t12 = Cummings PARSEC-based" << endl;
    cerr << "\t\t13 = Cummings MIST-based\n" << endl;

    cerr << "\t--wdModel" << endl;
    cerr << "\t\t0 = Wood" << endl;
    cerr << "\t\t1 = Montgomery (original)" << endl;
    cerr << "\t\t2 = Althaus" << endl;
    cerr << "\t\t3 = Renedo" << endl;
    cerr << "\t\t4 = Montgomery (2018)\n" << endl;

    cerr << "\t--bdModel" << endl;
    cerr << "\t\t0 = None" << endl;
    cerr << "\t\t1 = Baraffe\n" << endl;

    cerr << "\t--photFile" << endl;
    cerr << "\t\tThe absolute path to a photometry file for input\n" << endl;

    cerr << "\t--scatterFile" << endl;
    cerr << "\t\tThe absolute path to a scatter file for output. The scatter\n";
    cerr << "\t\tfile does not have to exist and WILL be overwritten without\n";
    cerr << "\t\twarning.\n" << endl;

    cerr << "\t--outputFileBase" << endl;
    cerr << "\t\tRun information is appended to this name\n" << endl;

    cerr << "\t--modelDirectory" << endl;
    cerr << "\t\tThe directory in which models are located\n" << endl;

    cerr << "\t--M_wd_up" << endl;
    cerr << "\t\tThe maximum mass for a WD-producing star\n" << endl;

    cerr << "\t--priorFe_H" << endl;
    cerr << "\t--sigmaFe_H" << endl;
    cerr << "\t--startingFe_H\n" << endl;

    cerr << "\t--priorDistMod" << endl;
    cerr << "\t--sigmaDistMod" << endl;
    cerr << "\t--startingDistMod\n" << endl;

    cerr << "\t--priorAv" << endl;
    cerr << "\t--sigmaAv" << endl;
    cerr << "\t--startingAv\n" << endl;

    cerr << "\t--priorY" << endl;
    cerr << "\t--sigmaY" << endl;
    cerr << "\t--startingY\n" << endl;

    cerr << "\t--priorCarbonicity" << endl;
    cerr << "\t--sigmaCarbonicity" << endl;
    cerr << "\t--startingCarbonicity\n" << endl;

    cerr << "\t--priorLogAge" << endl;
    cerr << "\t--sigmaLogAge" << endl;
    cerr << "\t--startingLogAge\n" << endl;

    cerr << "\t--minMag" << endl;
    cerr << "\t--maxMag\n" << endl;

    cerr << "\t--index" << endl;
    cerr << "\t\t0 being the first filter in the dataset\n" << endl;

    cerr << "\t--burnIter" << endl;
    cerr << "\t--maxIter\n" << endl;

    cerr << "\t--thin" << endl;
    cerr << "\t--nStars\n" << endl;

    cerr << "\t--percentBinary" << endl;
    cerr << "\t\tpercent binaries (drawn randomly)\n" << endl;

    cerr << "\t--percentDB" << endl;
    cerr << "\t\tpercent of WDs that have He atmospheres (drawn randomly)\n" << endl;

    cerr << "\t--nFieldStars\n" << endl;

    cerr << "\t--brightLimit" << endl;
    cerr << "\t\tapparant mags, can remove bright stars, e.g. RGB\n" << endl;

    cerr << "\t--faintLimit" << endl;
    cerr << "\t\tapparant mags, can remove faint stars, e.g. faint MS and WDs\n" << endl;

    cerr << "\t--relevantFilt" << endl;
    cerr << "\t\t0=bluest band available\n" << endl;

    cerr << "\t--limitS2N" << endl;
    cerr << "\t\tuse to remove objects with overly low signal-to-noise\n" << endl;

    cerr << "\t--deltaMass" << endl;
    cerr << "\t--deltaMassRatio" << endl;
    cerr << "\t\tSets the step sizes for mass and mass ratio in sampleMass\n" << endl;


    cerr << "\n9.3.0 flags" << endl;
    cerr <<   "===========" << endl;

    cerr << "\t--threads" << endl;
    cerr << "\t\tSpecify the number of local threads with which to run\n" << endl;

    cerr << "\t--bigStepBurnin" << endl;
    cerr << "\t\tRun the burnin only using the \"propClustBigSteps\" algorithm" << endl;


    cerr << "\n9.4.0 flags" << endl;
    cerr <<   "===========" << endl;

    cerr << "\t--noBinaries" << endl;
    cerr << "\t\tTurns off integration over secondary mass" << endl;


    cerr << "\n9.4.4 flags" << endl;
    cerr <<   "===========" << endl;

    cerr << "\t--noWDs" << endl;
    cerr << "\t\tKeeps simCluster from generating WDs" << endl;


    cerr << "\n9.5.0 flags" << endl;
    cerr <<   "===========" << endl;

    cerr << "\t--backend <int>" << endl;
    cerr << "\t\tSpecify the desired back end:" << endl;
    cerr << "\t\t\t0 = File" << endl;
    cerr << "\t\t\t1 = SQLite\n" << endl;

    cerr << "\t--priorParallax" << endl;
    cerr << "\t--sigmaParallax" << endl;
    cerr << "\t--startingParallax" << endl;
    cerr << "\t\tSpecifies distance values in parallax. Can not be combined\n";
    cerr << "\t\twith any DistMod flag. Must all be specified simultaneously." << endl;

    cerr << "\n9.6.0 flags" << endl;
    cerr << "\t--onlyWDs" << endl;
    cerr << "\t\tRestrict WD numerical integration to only occur above the AGB tip. Primarily useful when the MS model doesn't contain the filters you want in a WD-only run." << endl;

    cerr << "\n\t--allowInvalidModels" << endl;
    cerr << "\t\tAllows runs with MS models which would be invalid due to missing filters to proceed. WILL cause an unexplained crash if you actually use the MS models." << endl;


    cerr << "\n\t--wdAtmosphereModel" << endl;
    cerr << "\t\t0 = Bergeron (~2013)" << endl;
    cerr << "\t\t1 = Bergeron (2019)" << endl;
    cerr << "\t\t2 = Bergeron (2020)" << endl;

    cerr << "\n\t--eepInterpolationPower <int>" << endl;
    cerr << "\t\tAdjusts the number of interpolated steps between EEPs in photometry.\n" << endl;
    cerr << "\t\t       0 = Traditional EEP interpolation, 80 steps (default)" << endl;
    cerr << "\t\t-5 .. -1" << endl;
    cerr << "\t\t 1 ..  n = Configurable EEP interpolation, provides 2^(5 + n) steps." << endl;
    cerr << "\t\t\t   1 is on the order of the traditional value. Negative numbers offer a performance boost in exchange for precision and potentially accuracy." << endl;

    cerr << "\n\t--wdInterpolationPower <int>" << endl;
    cerr << "\t\tAs with eepInterpolationPower, but with different bounds for a different grid.\n" << endl;
    cerr << "\t\t        0 = Traditional WD interpolation, 8000 steps (default)" << endl;
    cerr << "\t\t-13 .. -1" << endl;
    cerr << "\t\t  1 ..  n = Configurable WD interpolation, provides 2^(12 + n) steps." << endl;
    cerr << "\t\t\t    1 is on the order of the traditional value. Negative numbers offer a performance boost in exchange for precision and potentially accuracy." << endl;

    cerr << "\n\t--includeBurnin" << endl;
    cerr << "\t\tReturns the pre-9.6 behavior to sampleMass and sampleWDMass where all lines from" << endl;
    cerr << "\t\tthe result file are considered (rather than just stage 3, post-burnin)." << endl;

    cerr << "\n\t--stopAfterBurnin" << endl;
    cerr << "\t\tStops the run after the Stage 2 burnin for use with startWithRes." << endl;

    cerr << "\n\t--startWithBurnin <filename>" << endl;
    cerr << "\t\tStarts the run with the Stage 2 burnin output from a standard run of BASE-9 or" << endl;
    cerr << "\t\ta run using the 'stopAfterBurnin' flag." << endl;

    cerr << "\n9.6.1 flags" << endl;
    cerr << "\n\t--exitAfterLoad" << endl;
    cerr << "\t\tDev flag for model and photometry debugging" << endl;

    cerr << "\n\t--ignoreLowSigma" << endl;
    cerr << "\t\tTurn off warning about sigmas in photometry less than 0.01" << endl;

    cerr << "\n\t--allowNegativeSigma" << endl;
    cerr << "\t\tAllow negative sigmas in photometry, indicating an ignored filter for that star" << endl;

    cerr << "\n\t--include LIST and --exclude LIST" << endl;
    cerr << "\t\tFilters the input photometry. LIST should be a comma-separated of IDs exactly as they are in the photometry (excluding leading spaces," << endl;
    cerr << "\t\tsee the output with `--veryVerbose --exitAfterLoad` flags if you are unsure)." << endl;

    cerr << "\n\t\tInvalid items in LIST (i.e., those not the photometry file) are ignored silently. The list of filtered IDs is printed if --verbose is set." << endl;

    cerr << "\n\t\tBehavior:" << endl;
    cerr << "\t\t\t--include excludes all stars but those listed." << endl;
    cerr << "\t\t\t--exclude includes all but those listed." << endl;
    cerr << "\t\t\t--include and --exclude together excludes all stars except those included in --include and not excluded in --exclude (i.e., --exclude trumps)." << endl;    

    cerr << "\n\t\tExample:" << endl;
    cerr << "\t\t\t$ ./bin/singlePopMcmc --photFile Hyades.UBV.phot --verbose --include \"vB043,HZ4,VR7\" --exclude \"HZ4\"" << endl;
    cerr << "\n\t\t\tOutput includes:" << endl;
    cerr << "\t\t\t\tLoaded stars: 'vB043', 'VR7'" << endl;

}

static void printVersion()
{
    cerr << "Base-" << Base9_VERSION << endl;
}
