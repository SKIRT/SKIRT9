/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SKIRTCOMMANDLINEHANDLER_HPP
#define SKIRTCOMMANDLINEHANDLER_HPP

#include "CommandLineArguments.hpp"
#include "ConsoleLog.hpp"

////////////////////////////////////////////////////////////////////

/**
This class processes the command line arguments for SKIRT and invokes the appropriate
high-level functions to perform the actions requested by the user. When invoked with invalid
command line arguments, it prints a brief help message. When invoked without any arguments, it
enters interactive mode, constructing a simulation from the user's responses and saving the
result in a ski file, without actually performing the simulation. Otherwise, it runs the
simulations in the ski files specified on the command line according to the following syntax:

\verbatim
 skirt [-t <threads>] [-s <simulations>] [-d]
       [-b] [-v] [-m] [-e]
       [-k] [-i <dirpath>] [-o <dirpath>]
       [-r] {<filepath>}*
\endverbatim

- The -t option specifies the number of parallel threads for each simulation. The default value
  is the number of logical cores on the computer running SKIRT.

- The -s option specifies the number of simulations to be executed in parallel. The default value is one.

- The -d option enables data parallelization mode for multiple processes.

- The -b option forces brief console logging, i.e. only success and error messages are shown rather than all progress
  messages. If there are multiple parallel simulations (see the -s option), the -b option is turned on automatically
  to avoid a plethora of randomly intermixing messages. If there is only one simulation at a time, the console shows
  all messages unless the -b option is present. In any case, the complete log output for each simulation is always
  written to a file in the output directory.

- The -v option enables verbose logging for simulations running with multiple processes, causing each process to
  create its own log file rather than relying on the root process to log all relevant information. The individual log
  files then also include messages bracketing each MPI operation, facilitating the diagnose of issues related to
  communication between the processes.

- The -m option causes information on current memory usage to be included in each log message.

- The -e option activates emulation mode, which can be used to estimate the amount of memory used by
  a given simulation without actually performing the simulation.

- The -k option causes the simulation input/output paths to be relative to the ski file being processed, rather than
  to the current directory. This is useful, for example, when processing multiple ski files organized in a nested
  directory hierarchy (see the -r option).

- The -i option specifies the absolute or relative path for simulation input files.

- The -o option specifies the absolute or relative path for simulation output files.

- The -r option causes recursive directory descent for all specified \<filepath\> arguments, in other words
  all directories inside the specified base paths are searched for the specified filename (or filename pattern).

In the simplest case, a \<filepath\> argument specifies the relative or absolute file path for a
single ski file, with or without the ".ski" filename extension. However the filename (\em not the base path)
may also contain ? and * wildcards forming a pattern to match multiple files. If the -r option
is present, all directories recursively nested within the base path are searched as well, using
the same filename pattern. If the filename contains wildcards or the -r option is present (in
other words, if the filepath may match multiple files) the ".ski" filename extension is not automatically added.
Furthermore, filepaths containing wildcards should be enclosed in quotes on the command
line to avoid expansion of the wildcard pattern by the shell.

For example, to process all "test" ski files inside the "geometry" directory hierarchy, one
could specify:

\verbatim
 skirt -s 4 -t 1 -r "/root-test-file-path/geometry/test*.ski"
\endverbatim
*/
class SkirtCommandLineHandler final
{
public:
    /** The constructor obtains the program's command line arguments and issues a welcome message
        to the console log. */
    SkirtCommandLineHandler();

    /** This function processes the command line arguments and invokes the appropriate high-level
        functions to perform the actions requested by the user. The function returns an appropriate
        application exit value. */
    int perform();

private:
    /** This function conducts an interactive session to construct a simulation and save the result
        in a ski file. The function returns an appropriate application exit value. */
    int doInteractive();

    /** This function scans the filepaths specified on the command line for ski files and performs
        the corresponding simulations according to the specified command line options. The function
        returns an appropriate application exit value. */
    int doBatch();

    /** This function exports a smile schema. This is an undocumented option. */
    int doSmileSchema();

    /** This function adds the ski filenames corresponding to the specified filepath to the
        internal list, after processing any wildcards and performing recursive descent if so
        requested by the -r option. If no appropriate filenames are found, the function logs a
        corresponding error message (and leaves the internal list unchanged). */
    void addSkiFilesFor(string filepath);

    /** This function adds the ski filenames corresponding to the specified name pattern inside the
        specified directory to the internal list. If so requested by the -r option, this function
        implements recursive descent by calling itself recursively for each subdirectory. */
    void addSkiFilesFor(string dirpath, string name);

    /** This function actually performs a single simulation constructed from the ski file at the
        specified index in the internal list. */
    void doSimulation(size_t index);

    /** This function logs a simulation construction error to an appropriate emergency log file
        with a name and location corresponding to the regular simulation log file. */
    void logErrorToFile(const vector<string>& message, string skipath);

    /** This function prints a brief help message to the console. */
    void printHelp();

private:
    // data members
    CommandLineArguments _args;
    ConsoleLog _console;
    string _producerInfo;
    string _hostUserInfo;
    vector<string> _skifiles;
    int _parallelSims{1};
    bool _hasError{false};
};

////////////////////////////////////////////////////////////////////

#endif
