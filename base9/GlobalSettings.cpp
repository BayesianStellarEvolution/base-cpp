#include <string>
#include <iostream>
#include <vector>

#include <cstdio>
#include <cstring>
#include <getopt.h>

#include "Base9Config.h"
#include "yaml-cpp/yaml.h"
#include "GlobalSettings.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::istringstream;
using std::string;
using std::vector;
using YAML::Node;
using YAML::LoadFile;

// Forward declaration
static void printVersion ();

// Helper function for parsing comma-separated lists
static void commaSeparated(string toParse, std::vector<std::string> &vec)
{
    std::stringstream ss(toParse);
    std::string lid;

    while (std::getline(ss, lid, ',')) {
        vec.push_back(lid);
    }
}


void GlobalSettings::loadSettings(int argc, char** argv, const string& defaultFile)
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


void GlobalSettings::fromYaml (const string& yamlFile)
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

    verbose = getOrRequest <int>(generalNode, "verbose");

    files.backend = static_cast<Backend>(getOrRequest <int>(filesNode, "backend"));
    files.phot = getOrRequest <string>(filesNode, "photFile");
    files.output = getOrRequest <string>(filesNode, "outputFileBase");
    files.scatter = getOrRequest <string>(filesNode, "scatterFile");
    files.models = getOrRequest <string>(filesNode, "modelDirectory");
}


void GlobalSettings::fromCLI (int argc, char **argv)
{
    char **t_argv = new char*[argc];

    for (int i = 0; i<argc; i++)
    {
        t_argv[i] = new char[strlen (argv[i]) + 1];

        strcpy (t_argv[i], argv[i]);
    }

    static struct option long_options[] = {
        // Model options
        {"msModelFile", required_argument, 0, 0xFE},
        {"ifmr", required_argument, 0, 0xFD},
        {"wdModel", required_argument, 0, 0xFC},
        {"M_wd_up", required_argument, 0, 0xFA},
        {"wdAtmosphereModel", required_argument, 0, 0xBD},

        // Cluster prior options
        {"priorFe_H", required_argument, 0, 0xF8},
        {"sigmaFe_H", required_argument, 0, 0xF7},
        {"priorDistMod", required_argument, 0, 0xF6},
        {"sigmaDistMod", required_argument, 0, 0xF5},
        {"priorAv", required_argument, 0, 0xF4},
        {"sigmaAv", required_argument, 0, 0xF3},
        {"priorY", required_argument, 0, 0xF2},
        {"sigmaY", required_argument, 0, 0xF1},
        {"priorLogAge", required_argument, 0, 0xF0},
        {"sigmaLogAge", required_argument, 0, 0xC9},
        {"priorCarbonicity", required_argument, 0, 0xCD},
        {"sigmaCarbonicity", required_argument, 0, 0xCC},
        {"priorParallax", required_argument, 0, 0xBF},
        {"sigmaParallax", required_argument, 0, 0xBE},

        // Cluster starting values
        {"startingFe_H", required_argument, 0, 0xC8},
        {"startingDistMod", required_argument, 0, 0xC7},
        {"startingAv", required_argument, 0, 0xC6},
        {"startingY", required_argument, 0, 0xC5},
        {"startingLogAge", required_argument, 0, 0xC4},
        {"startingCarbonicity", required_argument, 0, 0xC3},
        {"startingParallax", required_argument, 0, 0xC0},

        // Magnitude range
        {"minMag", required_argument, 0, 0xEF},
        {"maxMag", required_argument, 0, 0xEE},
        {"index", required_argument, 0, 0xED},

        // File options
        {"modelDirectory", required_argument, 0, 0xE5},
        {"photFile", required_argument, 0, 0xDF},
        {"scatterFile", required_argument, 0, 0xDE},
        {"outputFileBase", required_argument, 0, 0xDD},
        {"config", required_argument, 0, 0xDC},
        {"backend", required_argument, 0, 0xC2},

        // General options
        {"seed", required_argument, 0, 0xE0},
        {"help", no_argument, 0, 0xDB},
        {"version", no_argument, 0, 0xDA},

        // Flags
        {"verbose", no_argument, 0, 0xAF},
        {"veryVerbose", no_argument, 0, 0xB9},
        {"allowInvalidModels", no_argument, 0, 0x9E},
        {"details", no_argument, 0, 0xAA},
        {"ignoreLowSigma", no_argument, 0, 0x9C},
        {"allowNegativeSigma", no_argument, 0, 0x9B},
        {"includeBurnin", no_argument, 0, 0xBA},

        // Filter options
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

            case 0xC9:
                istringstream (string (optarg)) >> cluster.priorSigma.logAge;
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

            case 0xE5:
                files.models = optarg;
                break;

            case 0xE0:
                istringstream (string (optarg)) >> seed;
                break;

            case 0xDF:
                files.phot = optarg;
                break;

            case 0xDE:
                // --scatterFile is now app-specific (scatterCluster)
                break;

            case 0xDD:
                files.output = optarg;
                break;

            case 0xDC:
                files.config = optarg;
                break;

            case 0xDB:                  // --help
                GlobalSettings::printUsage ();
                exit (EXIT_SUCCESS);

            case 0xDA:
                printVersion ();
                exit (EXIT_SUCCESS);

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
                // --noBinaries is now app-specific
                break;

            case 0xAA:
                details = true;
                break;

            case 0x9F:
                // --onlyWDs is now app-specific
                break;

            case 0x9E:
                allowInvalidModels = true;
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
                // --eepInterpolationPower is now app-specific
                break;

            case 0xBB:
                // --wdInterpolationPower is now app-specific
                break;

            case 0xBA:
                includeBurnin = true;
                break;

            case 0x9A:
                commaSeparated(string(optarg), include);
                break;

            case 0x99:
                commaSeparated(string(optarg), exclude);
                break;

            case '?':
                // getopt_long already printed an error message.
                GlobalSettings::printUsage ();
                exit (EXIT_FAILURE);

            default:
                // Unknown option - ignore for app-specific options
                break;
        }
    }

    for (int i = 0; i<argc; i++)
    {
        delete[] t_argv[i];
    }
    delete[] t_argv;

    if (modBits == 0)  { return; }
    if (modBits == 7)  { modIsParallax = false; return; }
    if (modBits == 56) { modIsParallax = true;  return; }

    exitWith("You must specify starting, prior mean, and prior sigma values for either parallax or distMod.");
}

template<typename T>
T GlobalSettings::getDefault (Node & n, string f, T def)
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
T GlobalSettings::getOrDie (Node & n, string f)
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
T GlobalSettings::getOrRequest (Node & n, string f)
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

// Explicit template instantiations for common types
template int GlobalSettings::getDefault<int>(Node&, string, int);
template double GlobalSettings::getDefault<double>(Node&, string, double);
template bool GlobalSettings::getDefault<bool>(Node&, string, bool);
template string GlobalSettings::getDefault<string>(Node&, string, string);

template int GlobalSettings::getOrDie<int>(Node&, string);
template double GlobalSettings::getOrDie<double>(Node&, string);
template bool GlobalSettings::getOrDie<bool>(Node&, string);
template string GlobalSettings::getOrDie<string>(Node&, string);

template int GlobalSettings::getOrRequest<int>(Node&, string);
template double GlobalSettings::getOrRequest<double>(Node&, string);
template bool GlobalSettings::getOrRequest<bool>(Node&, string);
template string GlobalSettings::getOrRequest<string>(Node&, string);

Node GlobalSettings::getNode (Node & n, string f)
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

[[noreturn]] void GlobalSettings::exitWith (string s)
{
    cerr << s << endl;
    abort ();
}

void GlobalSettings::printUsage ()
{
    cerr << "\nBASE-9 Global Options\n";
    cerr << "=====================\n\n";

    cerr << "General Options:\n";
    cerr << "  --help                    Show this help message and exit\n";
    cerr << "  --version                 Show version information and exit\n";
    cerr << "  --config <file>           YAML configuration file (default: base9.yaml)\n";
    cerr << "  --seed <int>              Random number generator seed\n";
    cerr << "  --verbose                 Enable verbose output\n";
    cerr << "  --veryVerbose             Enable very verbose output\n";
    cerr << "\n";

    cerr << "File Options:\n";
    cerr << "  --photFile <path>         Input photometry file path\n";
    cerr << "  --outputFileBase <path>   Output file base name\n";
    cerr << "  --modelDirectory <path>   Directory containing stellar models\n";
    cerr << "  --backend <int>           Output backend: 0=File, 1=SQLite\n";
    cerr << "\n";

    cerr << "Main Sequence Model:\n";
    cerr << "  --msModelFile <file>      Model file relative to modelDirectory\n";
    cerr << "\n";

    cerr << "White Dwarf Model:\n";
    cerr << "  --ifmr <int>              Initial-final mass relation:\n";
    cerr << "                              0=Weidemann, 1=Williams, 2=Salaris linear,\n";
    cerr << "                              3=Salaris piecewise linear, 4-11=tunable,\n";
    cerr << "                              12=Cummings PARSEC, 13=Cummings MIST\n";
    cerr << "  --wdModel <int>           WD cooling model:\n";
    cerr << "                              0=Wood, 1=Montgomery (original),\n";
    cerr << "                              2=Althaus, 3=Renedo, 4=Montgomery (2018)\n";
    cerr << "  --wdAtmosphereModel <int> WD atmosphere model:\n";
    cerr << "                              0=Bergeron (~2013), 1=Bergeron (2019),\n";
    cerr << "                              2=Bergeron (2020), 3=Bergeron (2026 HST)\n";
    cerr << "  --M_wd_up <float>         Maximum mass for WD-producing star\n";
    cerr << "\n";

    cerr << "Cluster Parameters - Starting Values:\n";
    cerr << "  --startingFe_H <float>         Starting [Fe/H]\n";
    cerr << "  --startingDistMod <float>      Starting distance modulus\n";
    cerr << "  --startingParallax <float>     Starting parallax (mas)\n";
    cerr << "  --startingAv <float>           Starting extinction Av\n";
    cerr << "  --startingY <float>            Starting helium abundance Y\n";
    cerr << "  --startingLogAge <float>       Starting log(age) in years\n";
    cerr << "  --startingCarbonicity <float>  Starting carbonicity\n";
    cerr << "\n";

    cerr << "Cluster Parameters - Prior Means:\n";
    cerr << "  --priorFe_H <float>            Prior mean [Fe/H]\n";
    cerr << "  --priorDistMod <float>         Prior mean distance modulus\n";
    cerr << "  --priorParallax <float>        Prior mean parallax (mas)\n";
    cerr << "  --priorAv <float>              Prior mean extinction Av\n";
    cerr << "  --priorY <float>               Prior mean helium abundance Y\n";
    cerr << "  --priorLogAge <float>          Prior mean log(age) in years\n";
    cerr << "  --priorCarbonicity <float>     Prior mean carbonicity\n";
    cerr << "\n";

    cerr << "Cluster Parameters - Prior Sigmas:\n";
    cerr << "  --sigmaFe_H <float>            Prior sigma [Fe/H]\n";
    cerr << "  --sigmaDistMod <float>         Prior sigma distance modulus\n";
    cerr << "  --sigmaParallax <float>        Prior sigma parallax (mas)\n";
    cerr << "  --sigmaAv <float>              Prior sigma extinction Av\n";
    cerr << "  --sigmaY <float>               Prior sigma helium abundance Y\n";
    cerr << "  --sigmaLogAge <float>          Prior sigma log(age)\n";
    cerr << "  --sigmaCarbonicity <float>     Prior sigma carbonicity\n";
    cerr << "\n";

    cerr << "NOTE: Use either distMod OR parallax options, not both.\n";
    cerr << "      All three (starting, prior, sigma) must be specified for the chosen type.\n";
    cerr << "\n";

    cerr << "Magnitude Range:\n";
    cerr << "  --minMag <float>          Minimum magnitude\n";
    cerr << "  --maxMag <float>          Maximum magnitude\n";
    cerr << "  --index <int>             Filter index (0 = first filter in dataset)\n";
    cerr << "\n";


    cerr << "Star Filtering:\n";
    cerr << "  --include <list>          Include only these stars (comma-separated IDs)\n";
    cerr << "  --exclude <list>          Exclude these stars (comma-separated IDs)\n";
    cerr << "                            When both used, --exclude takes precedence\n";
    cerr << "\n";

    cerr << "Behavior Flags:\n";
    cerr << "  --allowInvalidModels      Allow models with missing filters\n";
    cerr << "  --includeBurnin           Include burn-in samples in output\n";
    cerr << "\n";

    cerr << "Debugging:\n";
    cerr << "  --details                 Show additional run details\n";
    cerr << "  --ignoreLowSigma          Suppress warnings for sigma < 0.01\n";
    cerr << "  --allowNegativeSigma      Allow negative sigmas (disables filter)\n";
    cerr << "\n";
}

void GlobalSettings::printVerboseUsage ()
{
    cerr << "\nDetailed Option Reference\n";
    cerr << "=========================\n\n";

    cerr << "Model Options:\n";
    cerr << "  --msModelFile <file>\n";
    cerr << "      Specify a model file relative to --modelDirectory\n\n";

    cerr << "  --ifmr <int>              Initial-final mass relation:\n";
    cerr << "      0 = Weidemann\n";
    cerr << "      1 = Williams\n";
    cerr << "      2 = Salaris linear\n";
    cerr << "      3 = Salaris piecewise linear\n";
    cerr << "      4-11 = Tunable\n";
    cerr << "      12 = Cummings PARSEC-based\n";
    cerr << "      13 = Cummings MIST-based\n\n";

    cerr << "  --wdModel <int>           WD cooling model:\n";
    cerr << "      0 = Wood\n";
    cerr << "      1 = Montgomery (original)\n";
    cerr << "      2 = Althaus\n";
    cerr << "      3 = Renedo\n";
    cerr << "      4 = Montgomery (2018)\n\n";

    cerr << "  --wdAtmosphereModel <int> WD atmosphere model:\n";
    cerr << "      0 = Bergeron (~2013)\n";
    cerr << "      1 = Bergeron (2019)\n";
    cerr << "      2 = Bergeron (2020)\n";
    cerr << "      3 = Bergeron (2026 [HST-only])\n\n";

    cerr << "File Options:\n";
    cerr << "  --scatterFile <path>\n";
    cerr << "      The absolute path to a scatter file for output. The scatter\n";
    cerr << "      file does not have to exist and WILL be overwritten without\n";
    cerr << "      warning.\n\n";

    cerr << "Cluster Parameters:\n";
    cerr << "  NOTE: Use either distMod OR parallax options, not both.\n";
    cerr << "        All three (starting, prior, sigma) must be specified for the chosen type.\n\n";

    cerr << "  --index <int>\n";
    cerr << "      0 being the first filter in the dataset\n\n";

    cerr << "Interpolation Power:\n";
    cerr << "  --eepInterpolationPower <int>\n";
    cerr << "      Adjusts the number of interpolated steps between EEPs in the isochrone.\n";
    cerr << "      The isochrone contains ~100 EEPs; this multiplies that by the step count.\n\n";
    cerr << "             0 = Traditional EEP interpolation, 80 steps per EEP interval (default)\n";
    cerr << "      -5 .. -1 = Reduced interpolation; -5 gives 1 step (un-interpolated model)\n";
    cerr << "       1 ..  n = Configurable EEP interpolation, provides 2^(5+n) steps per EEP interval.\n";
    cerr << "                 1 (64 steps) is on the order of the traditional value.\n\n";

    cerr << "  --wdInterpolationPower <int>\n";
    cerr << "      Adjusts the number of WD mass grid points between MIN_MASS1 and M_wd_up.\n\n";
    cerr << "             0 = Traditional WD interpolation, 8000 steps (default)\n";
    cerr << "      -7 .. -1 = Reduced interpolation; -7 gives 32 steps\n";
    cerr << "       1 ..  n = Configurable WD interpolation, provides 2^(12+n) steps.\n";
    cerr << "                 1 (8192 steps) is on the order of the traditional value.\n\n";

    cerr << "Star Filtering:\n";
    cerr << "  --include LIST and --exclude LIST\n";
    cerr << "      Filters the input photometry. LIST should be a comma-separated list of IDs\n";
    cerr << "      exactly as they are in the photometry (excluding leading spaces; see output\n";
    cerr << "      with --veryVerbose --exitAfterLoad flags if unsure).\n\n";

    cerr << "      Invalid items in LIST (i.e., those not in the photometry file) are ignored\n";
    cerr << "      silently. The list of filtered IDs is printed if --verbose is set.\n\n";

    cerr << "      Behavior:\n";
    cerr << "        --include excludes all stars but those listed.\n";
    cerr << "        --exclude includes all but those listed.\n";
    cerr << "        --include and --exclude together excludes all stars except those\n";
    cerr << "          included in --include and not excluded in --exclude.\n\n";

    cerr << "      Example:\n";
    cerr << "        $ ./bin/singlePopMcmc --photFile Hyades.UBV.phot --verbose \\\n";
    cerr << "            --include \"vB043,HZ4,VR7\" --exclude \"HZ4\"\n";
    cerr << "        Output includes: Loaded stars: 'vB043', 'VR7'\n\n";

    cerr << "Behavior Flags:\n";
    cerr << "  --onlyWDs\n";
    cerr << "      Restrict WD numerical integration to only occur above the AGB tip.\n";
    cerr << "      Primarily useful when the MS model doesn't contain the filters you want\n";
    cerr << "      in a WD-only run.\n\n";

    cerr << "  --allowInvalidModels\n";
    cerr << "      Allows runs with MS models which would be invalid due to missing filters\n";
    cerr << "      to proceed. WILL cause an unexplained crash if you actually use the MS models.\n\n";

    cerr << "  --includeBurnin\n";
    cerr << "      Returns the pre-9.6 behavior to sampleMass and sampleWDMass where all lines\n";
    cerr << "      from the result file are considered (rather than just stage 3, post-burnin).\n\n";

    cerr << "  --ignoreLowSigma\n";
    cerr << "      Turn off warning about sigmas in photometry less than 0.01\n\n";

    cerr << "  --allowNegativeSigma\n";
    cerr << "      Allow negative sigmas in photometry, indicating an ignored filter for that star\n\n";
}

static void printVersion()
{
    cerr << "Base-" << Base9_VERSION << endl;
}
