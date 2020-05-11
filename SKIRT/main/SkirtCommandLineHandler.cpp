/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SkirtCommandLineHandler.hpp"
#include "BuildInfo.hpp"
#include "Console.hpp"
#include "ConsoleHierarchyCreator.hpp"
#include "FatalError.hpp"
#include "FileLog.hpp"
#include "FilePaths.hpp"
#include "MonteCarloSimulation.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "ProcessManager.hpp"
#include "SchemaDef.hpp"
#include "SimulationItemRegistry.hpp"
#include "StopWatch.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "TimeLogger.hpp"
#include "XmlHierarchyCreator.hpp"
#include "XmlHierarchyWriter.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // the allowed options list, in the format consumed by the CommandLineArguments constructor
    static const char* allowedOptions = "-t* -s* -d -b -v -m -e -k -i* -o* -r -x";
}

////////////////////////////////////////////////////////////////////

SkirtCommandLineHandler::SkirtCommandLineHandler() : _args(System::arguments(), allowedOptions)
{
    // issue welcome message
    _producerInfo =
        "SKIRT " + BuildInfo::projectVersion() + " (" + BuildInfo::codeVersion() + " " + BuildInfo::timestamp() + ")";
    _hostUserInfo = "Running on " + System::hostname() + " for " + System::username();
    _console.info("Welcome to " + _producerInfo);
    _console.info(_hostUserInfo);
}

////////////////////////////////////////////////////////////////////

int SkirtCommandLineHandler::perform()
{
    // catch and properly report any exceptions
    try
    {
        // if there are no arguments at all --> interactive mode
        // if there is at least one file path argument --> batch mode
        // if the -x option is present --> export smile schema (undocumented option)
        // otherwise --> error
        if (_args.isValid() && !_args.hasOptions() && !_args.hasFilepaths()) return doInteractive();
        if (_args.hasFilepaths()) return doBatch();
        if (_args.isPresent("-x")) return doSmileSchema();
        _console.error("Invalid command line arguments");
        printHelp();
        return EXIT_FAILURE;
    }
    catch (FatalError& error)
    {
        for (string line : error.message()) _console.error(line);
    }
    catch (const std::exception& except)
    {
        _console.error("Standard Library Exception: " + string(except.what()));
    }
    ProcessManager::abort(EXIT_FAILURE);
    return EXIT_FAILURE;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // issue a log warning message if the Core resource pack has not been installed
    // or if the installed version of any resource pack does not match the expected version
    // has an unexpected version
    void reportResourceIssues(Log* log)
    {
        bool haveIssue = false;

        if (FilePaths::installedPackVersion("Core") <= 0)
        {
            log->warning("Core resource files have not been installed");
            haveIssue = true;
        }

        if (!haveIssue)
            for (string packname : FilePaths::expectedPacks())
            {
                int expected = FilePaths::expectedPackVersion(packname);
                int installed = FilePaths::installedPackVersion(packname);

                if (expected && installed && expected != installed)
                {
                    log->warning("Version number does not match for resource pack " + packname + ": expected "
                                 + std::to_string(expected) + ", installed " + std::to_string(installed));
                    haveIssue = true;
                }
            }

        if (haveIssue)
        {
            log->warning("  - run  './downloadResources.sh' from the ~/SKIRT/git directory");
            log->warning("  - refer to installation guide on www.skirt.ugent.be for more info");
        }
    }

    // issue a log info message reporting the peak memory usage so far
    void reportPeakMemory(Log* log)
    {
        size_t avail = System::availableMemory();
        size_t peak = System::peakMemoryUsage();
        double peakPerCent = 100. * static_cast<double>(peak) / static_cast<double>(avail);

        log->info("Available memory: " + StringUtils::toMemSizeString(avail) + " -- Peak memory usage: "
                  + StringUtils::toMemSizeString(peak) + " (" + StringUtils::toString(peakPerCent, 'f', 1) + "%)");
    }
}

////////////////////////////////////////////////////////////////////

int SkirtCommandLineHandler::doInteractive()
{
    if (ProcessManager::isMultiProc()) throw FATALERROR("Interactive mode cannot be run with multiple processes");

    // alert the user about problems with the installed resource packs
    reportResourceIssues(&_console);

    // ask for the name of the ski file in which to save the result
    Console::info("Interactively constructing a simulation...");
    string filename;
    while (true)
    {
        // get a file name, adding the .ski extension if needed
        filename = Console::promptForString("Enter the name of the ski file to be created", false, string());
        filename = StringUtils::addExtension(filename, "ski");

        // verify that the file does not yet exist
        // (we test whether the file can be opened, which is the best we can do in standard C++14)
        if (System::ifstream(filename))
            Console::error("This file already exists; enter another name");
        else
            break;
    }

    // interactively construct the simulation item hierarchy
    auto schema = SimulationItemRegistry::getSchemaDef();
    auto simulation = ConsoleHierarchyCreator::create(schema);  // unique pointer

    // create the ski file reflecting this simulation
    XmlHierarchyWriter::write(simulation.get(), schema, filename, _producerInfo);

    // report success
    Console::info("Successfully created ski file '" + filename + "'.");
    Console::info("To run the simulation use the command: skirt " + filename.substr(0, filename.length() - 4));
    return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////

int SkirtCommandLineHandler::doBatch()
{
    // build a list of filenames for existing ski files
    _skifiles.clear();
    _hasError = false;
    for (string filepath : _args.filepaths()) addSkiFilesFor(filepath);

    // exit if there were any problems with the file paths
    if (_hasError)
    {
        if (!_args.isPresent("-b")) printHelp();
        return EXIT_FAILURE;
    }

    // if there is only one ski file, simply perform the single simulation
    size_t numSkiFiles = _skifiles.size();
    if (numSkiFiles == 1)
    {
        _parallelSims = 1;
        doSimulation(0);
    }
    else
    {
        // determine the number of parallel simulations
        _parallelSims = max(_args.intValue("-s"), 1);

        // handle the serial case separately to avoid using MPI nested within a Parallel instance
        if (_parallelSims == 1)
        {
            // perform a simulation for each ski file
            TimeLogger logger(&_console, "a set of " + std::to_string(numSkiFiles) + " simulations");
            for (size_t i = 0; i != numSkiFiles; ++i) doSimulation(i);
        }
        else
        {
            // prevent multiple simulations to be launched in parallel while MPI parallelization is used
            if (ProcessManager::isMultiProc())
                throw FATALERROR("Cannot run multiple simulations in parallel when there are multiple MPI processes");

            // perform a simulation for each ski file
            TimeLogger logger(&_console, "a set of " + std::to_string(numSkiFiles) + " simulations, "
                                             + std::to_string(_parallelSims) + " in parallel");
            ParallelFactory factory;
            factory.setMaxThreadCount(_parallelSims);
            factory.parallelRootOnly()->call(numSkiFiles, [this](size_t first, size_t size) {
                for (size_t i = 0; i != size; ++i) doSimulation(first + i);
            });
        }
    }

    // report memory statistics for the complete run
    reportPeakMemory(&_console);

    // report stopwatch results, if any
    for (string line : StopWatch::report()) _console.warning(line);
    return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////

int SkirtCommandLineHandler::doSmileSchema()
{
    auto schema = SimulationItemRegistry::getSchemaDef();
    schema->save("skirt.smile", _producerInfo);
    _console.info("Successfully created SMILE schema file 'skirt.smile'.");
    return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////

void SkirtCommandLineHandler::addSkiFilesFor(string filepath)
{
    string name = StringUtils::filename(filepath);

    // no recursion and no wildcards -> expect a single result
    if (!_args.isPresent("-r") && !StringUtils::contains(name, "?") && !StringUtils::contains(name, "*"))
    {
        // if the file does not exist as specified, try adding the .ski extension
        if (!System::ifstream(filepath)) filepath = StringUtils::addExtension(filepath, "ski");
        if (!System::ifstream(filepath))
        {
            _console.error("This ski file does not exist: " + filepath);
            _hasError = true;
        }
        else
        {
            _skifiles.push_back(filepath);
        }
    }

    // recursion and/or wildcards -> multiple results possible
    else
    {
        // find matching files in this directory, possibly recursively (depending on -r option);
        // in this case do not automatically add the .ski extension; it leads to trouble with patterns
        // like "test*" which are automatically expanded by the shell before invoking the application
        auto oldSize = _skifiles.size();
        addSkiFilesFor(StringUtils::dirPath(filepath), name);
        auto newSize = _skifiles.size();
        if (newSize == oldSize)
        {
            _hasError = true;
            _console.error("No ski file matches the pattern: " + filepath);
        }
    }
}

////////////////////////////////////////////////////////////////////

void SkirtCommandLineHandler::addSkiFilesFor(string dirpath, string name)
{
    // add matching files at the current directory level
    for (string candidate : System::filesInDirectory(dirpath))
    {
        if (StringUtils::matches(candidate, name)) _skifiles.push_back(StringUtils::joinPaths(dirpath, candidate));
    }

    // if recursion is requested, call ourselves for all directories at this level
    if (_args.isPresent("-r"))
    {
        for (string subdir : System::dirsInDirectory(dirpath))
        {
            addSkiFilesFor(StringUtils::joinPaths(dirpath, subdir), name);
        }
    }
}

////////////////////////////////////////////////////////////////////

void SkirtCommandLineHandler::doSimulation(size_t index)
{
    if (_skifiles.size() > 1)
        _console.warning("Performing simulation #" + std::to_string(index + 1) + " of "
                         + std::to_string(_skifiles.size()));
    string skipath = _skifiles[index];
    _console.info("Constructing a simulation from ski file '" + skipath + "'...");

    // flag becomes true as soon as the simulation log file is available and used for reporting errors
    bool running = false;

    // construct and run the simulation; catch and rethrow exceptions so they are also logged to file
    try
    {
        // construct the simulation hierarchy from the ski file
        auto schema = SimulationItemRegistry::getSchemaDef();
        auto topitem = XmlHierarchyCreator::readFile(schema, skipath);  // unique pointer to Item
        auto simulation = dynamic_cast<MonteCarloSimulation*>(topitem.get());

        // set up simulation attributes that are not loaded from the ski file:
        //  - the paths for input and output files
        simulation->filePaths()->setOutputPrefix(StringUtils::filenameBase(skipath));
        string base = _args.isPresent("-k") ? StringUtils::dirPath(skipath) : "";
        string inpath = _args.value("-i");
        string outpath = _args.value("-o");
        if (!StringUtils::isAbsolutePath(inpath)) inpath = StringUtils::joinPaths(base, inpath);
        if (!StringUtils::isAbsolutePath(outpath)) outpath = StringUtils::joinPaths(base, outpath);
        simulation->filePaths()->setInputPath(inpath);
        simulation->filePaths()->setOutputPath(outpath);

        //  - the number of parallel threads
        if (_args.intValue("-t") > 0) simulation->parallelFactory()->setMaxThreadCount(_args.intValue("-t"));

        //  - the activation of data parallelization
        if (_args.isPresent("-d") && ProcessManager::isMultiProc())
        {
            throw FATALERROR("Data parallelization (-d option) is not supported at this time");
        }

        //  - the logging mechanisms
        FileLog* log = new FileLog();
        simulation->log()->setLinkedLog(log);
        simulation->log()->setVerbose(_args.isPresent("-v"));
        simulation->log()->setMemoryLogging(_args.isPresent("-m"));
        if (_parallelSims > 1 || _args.isPresent("-b")) simulation->log()->setLowestLevel(Log::Level::Success);

        // output a ski file reflecting this simulation for later reference
        if (ProcessManager::isRoot())
        {
            auto schema = SimulationItemRegistry::getSchemaDef();
            XmlHierarchyWriter::write(simulation, schema, simulation->filePaths()->output("parameters.xml"),
                                      _producerInfo);
        }

        // put the simulation in emulation mode if requested
        if (_args.isPresent("-e"))
        {
            simulation->log()->setLowestLevel(Log::Level::Error);
            simulation->config()->setEmulationMode();
        }

        // issue welcome message to the simulation log file
        log->setup();
        log->info(_producerInfo);
        log->info(_hostUserInfo);

        // log a warning about problems with the installed resource packs
        reportResourceIssues(simulation->log());

        // run the simulation and catch and properly report any exceptions to the similation log file
        try
        {
            running = true;
            simulation->setupAndRun();
        }
        catch (FatalError& error)
        {
            for (string line : error.message()) log->error(line);
            throw error;
        }
        catch (const std::exception& except)
        {
            log->error("Standard Library Exception: " + string(except.what()));
            throw except;
        }

        // if this is the only or first simulation in the run, report memory statistics in the simulation's log file
        if (_parallelSims == 1 && index == 0) reportPeakMemory(_args.isPresent("-v") ? simulation->log() : log);
    }
    catch (FatalError& error)
    {
        if (!running) logErrorToFile(error.message(), skipath);
        throw error;
    }
    catch (const std::exception& except)
    {
        if (!running) logErrorToFile(vector<string>({"Standard Library Exception: " + string(except.what())}), skipath);
        throw except;
    }
}

////////////////////////////////////////////////////////////////////

void SkirtCommandLineHandler::logErrorToFile(const vector<string>& message, string skipath)
{
    // construct the log file path
    string prefix = StringUtils::filenameBase(skipath);
    string base = _args.isPresent("-k") ? StringUtils::dirPath(skipath) : "";
    string outpath = _args.value("-o");
    if (!StringUtils::isAbsolutePath(outpath)) outpath = StringUtils::joinPaths(base, outpath);
    string logpath = StringUtils::joinPaths(outpath, prefix + "_log.txt");

    // create the log file
    std::ofstream logfile = System::ofstream(logpath);
    string stamp1 = System::timestamp() + "   ";
    logfile << stamp1 << _producerInfo << std::endl;
    logfile << stamp1 << _hostUserInfo << std::endl;
    string stamp2 = System::timestamp() + " * *** Error: ";
    logfile << stamp2 << "A fatal error occurred while constructing a simulation" << std::endl;
    logfile << stamp2 << "From ski file " + System::canonicalPath(skipath) << std::endl;
    for (string line : message) logfile << stamp2 << line << std::endl;
}

////////////////////////////////////////////////////////////////////

void SkirtCommandLineHandler::printHelp()
{
    if (!ProcessManager::isRoot()) return;

    _console.warning("");
    _console.warning("To create a new ski file interactively:    skirt");
    _console.warning("To run a simulation with default options:  skirt <ski-filename>");
    _console.warning("");
    _console.warning("  skirt [-t <threads>] [-s <simulations>] [-d]");
    _console.warning("        [-b] [-v] [-m] [-e]");
    _console.warning("        [-k] [-i <dirpath>] [-o <dirpath>]");
    _console.warning("        [-r] {<filepath>}*");
    _console.warning("");
    _console.warning("  -t <threads> : the number of parallel threads for each simulation");
    _console.warning("  -s <simulations> : the number of parallel simulations per process");
    _console.warning("  -d : enable data parallelization mode for multiple processes");
    _console.warning("  -b : force brief console logging");
    _console.warning("  -v : force verbose logging for multiple processes");
    _console.warning("  -m : state the amount of used memory at the start of each log message");
    _console.warning("  -e : run the simulation in emulation mode to get an estimate of the memory consumption");
    _console.warning("  -k : make the input/output paths relative to the ski file being processed");
    _console.warning("  -i <dirpath> : the relative or absolute path for simulation input files");
    _console.warning("  -o <dirpath> : the relative or absolute path for simulation output files");
    _console.warning("  -r : cause recursive directory descent for all specified ski file paths");
    _console.warning("  <filepath> : the relative or absolute file path for a ski file");
    _console.warning("               (the filename may contain ? and * wildcards)");
    _console.warning("");
}

//////////////////////////////////////////////////////////////////////
