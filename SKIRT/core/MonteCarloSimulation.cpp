/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MonteCarloSimulation.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "TimeLogger.hpp"

// included for testing purposes
#include "FilePaths.hpp"

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setupSelfBefore()
{
    Simulation::setupSelfBefore();
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setEmulationMode()
{
    _emulationMode = true;
    _numPackages = 0;
}

////////////////////////////////////////////////////////////////////

bool MonteCarloSimulation::emulationMode()
{
    return _emulationMode;
}

////////////////////////////////////////////////////////////////////

int MonteCarloSimulation::dimension() const
{
    return 0;
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runSelf()
{
    TimeLogger logger(log(), "the test phase");

    log()->info("Located resource "+ FilePaths::resource("README.txt"));
    log()->info("Located resource "+ FilePaths::resource("SunSED.stab"));

    auto map = System::acquireMemoryMap(FilePaths::resource("SunSED.stab"));
    auto map2 = System::acquireMemoryMap(FilePaths::resource("SunSED.stab"));
    if (map != map2) log()->error("Mappings differ !");

    auto start = static_cast<const char*>(map.first);
    auto end = static_cast<const char*>(map.first) + map.second;

    if (map.first)
    {
        log()->warning(string(start, 7));
        log()->warning(string(end-8, 7));
    }
}

////////////////////////////////////////////////////////////////////
