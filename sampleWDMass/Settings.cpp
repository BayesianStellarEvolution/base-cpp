#include <string>
#include <iostream>

#include <cstring>
#include <getopt.h>

#include "Settings.hpp"

using std::cerr;
using std::endl;
using std::string;

// Application-specific options (file scope for access from loadSettings)
static struct option long_options[] = {
    {"onlyWDs", no_argument, 0, 0x9F},
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
    }
    else
    {
        global.fromYaml(defaultFile);
    }

    global.fromCLI(argc, argv, long_options);
    fromCLI(argc, argv);
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
            case 0x9F:
                onlyWDs = true;
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
    cerr << "\nsampleWDMass - Sample White Dwarf Masses\n";
    cerr << "=========================================\n\n";

    cerr << "Samples white dwarf masses from posterior distributions\n";
    cerr << "given cluster parameters and observed photometry.\n\n";

    cerr << "Application-Specific Options:\n";
    cerr << "  --onlyWDs                 Integrate only above AGB tip (WD-only runs)\n";
    cerr << "\n";

}
