/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Simulation.hpp"
#include "FatalError.hpp"
#include "ProcessManager.hpp"
#include "TimeLogger.hpp"

////////////////////////////////////////////////////////////////////

void Simulation::setup()
{
    if (setupStarted()) return;

    _log->setup();
    TimeLogger logger(_log, "setup");
    SimulationItem::setup();

    find<Log>()->info("Waiting for other processes to finish setup...");
    ProcessManager::wait();
}

////////////////////////////////////////////////////////////////////

void Simulation::run()
{
    // verify setup
    if (!setupStarted()) throw FATALERROR("Simulation has not been setup before being run");

    // make each process use different random seeds
    if (ProcessManager::isMultiProc()) _random->randomize();

    TimeLogger logger(_log, "the simulation run");
    runSelf();
}

////////////////////////////////////////////////////////////////////

void Simulation::setupAndRun()
{
    string processInfo;

    _factory->setup();
    auto threads = _factory->maxThreadCount();
    if (threads > 1) processInfo += " using " + std::to_string(threads) + " threads";
    else processInfo += " using a single thread";

    auto procs = ProcessManager::size();
    if (procs > 1) processInfo += " for each of " + std::to_string(procs) + " processes";
    else processInfo += " and a single process";

    _log->setup();
    TimeLogger logger(_log, "simulation " + _paths->outputPrefix() + processInfo);

    setup();
    run();
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
