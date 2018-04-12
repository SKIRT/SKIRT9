/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Simulation.hpp"
#include "FatalError.hpp"
#include "TimeLogger.hpp"

////////////////////////////////////////////////////////////////////

void Simulation::setup()
{
    if (setupStarted()) return;

    _log->setup();
    TimeLogger logger(_log, "setup");
    SimulationItem::setup();

    // Wait for the other processes to reach this point
    _communicator->wait("the setup");
}

////////////////////////////////////////////////////////////////////

void Simulation::run()
{
    // verify setup
    if (!setupStarted()) throw FATALERROR("Simulation has not been setup before being run");

    if (_communicator->isMultiProc()) _random->randomize();

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

    _communicator->setup();
    auto procs = _communicator->size();
    if (procs > 1)
    {
        processInfo += " for each of " + std::to_string(procs) + " processes";
        if (_communicator->dataParallel()) processInfo += " in data parallelization mode";
        else processInfo += " in task parallelization mode";
    }
    else
    {
        processInfo += " and a single process";
    }

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

PeerToPeerCommunicator* Simulation::communicator() const
{
    return _communicator;
}

////////////////////////////////////////////////////////////////////
