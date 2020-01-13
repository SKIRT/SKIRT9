/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Simulation.hpp"
#include "ProcessManager.hpp"
#include "TimeLogger.hpp"

////////////////////////////////////////////////////////////////////

void Simulation::setupAndRun()
{
    // get information on the number of executions threads and processes
    string processInfo;

    _factory->setup();
    auto threads = _factory->maxThreadCount();
    if (threads > 1)
        processInfo += " using " + std::to_string(threads) + " threads";
    else
        processInfo += " using a single thread";

    auto procs = ProcessManager::size();
    if (procs > 1)
        processInfo += " for each of " + std::to_string(procs) + " processes";
    else
        processInfo += " and a single process";

    // log a start/finish message (ensure that the logger is initialized)
    _log->setup();
    TimeLogger logger(_log, "simulation " + _paths->outputPrefix() + processInfo);

    // setup and run the simulation
    setupSimulation();
    runSimulation();
}

////////////////////////////////////////////////////////////////////

FilePaths* Simulation::filePaths() const
{
    return _paths;
}

////////////////////////////////////////////////////////////////////

Log* Simulation::log() const
{
    return _log;
}

////////////////////////////////////////////////////////////////////

ParallelFactory* Simulation::parallelFactory() const
{
    return _factory;
}

////////////////////////////////////////////////////////////////////
