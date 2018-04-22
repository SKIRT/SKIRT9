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
#include "StaticPaths.hpp"
#include "StoredTable.hpp"

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

    StoredTable<1> table;
    table.open(this, "SunSED", "lambda(m)", "Llambda(W/m)");

    log()->warning("0: " + StringUtils::toString(table.atIndices(0)));
    log()->warning("1: " + StringUtils::toString(table.atIndices(10)));
    log()->warning("2: " + StringUtils::toString(table.atIndices(700)));
}

////////////////////////////////////////////////////////////////////
