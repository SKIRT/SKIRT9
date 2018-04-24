/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MonteCarloSimulation.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "TimeLogger.hpp"

// included for testing purposes
#include "StoredTable.hpp"
#include "TextOutFile.hpp"

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

    {
        StoredTable<1> table;
        table.open(this, "SunSED", "lambda(m)", "Llambda(W/m)");

        log()->warning("0: " + StringUtils::toString(table[0.55e-6], 'e', 9));
        log()->warning("1: " + StringUtils::toString(table[1e-6], 'e', 9));
        log()->warning("2: " + StringUtils::toString(table[100e-6], 'e', 9));
    }
    {
        StoredTable<2> table;
        table.open(this, "DustEM_Gra_OpticalProps", "lambda(m),a(m)", "Qabs(1)");

        log()->warning("3: " + StringUtils::toString(table(0.55e-6, 0.01e-6), 'e', 9));
    }
    {
        StoredTable<1> table;
        table.open(this, "MeanIvezicBenchmarkOpticalProps", "lambda(m)", "sigmaabs(m2/H)");

        Array xv, Yv;
        log()->warning("4: " + StringUtils::toString(table.cdf(xv, Yv, 100, 0.1e-6, 1e-6), 'e', 9));
        log()->warning("5: " + StringUtils::toString(table.cdf(xv, Yv, 1000, 1e-6, 10e-6), 'e', 9));
    }
    {
        StoredTable<3> table;
        table.open(this, "BruzualCharlotSEDFamily_Chabrier_hr", "lambda(m),Z(1),t(yr)", "Llambda(W/m)");

        Array xv, Yv;
        log()->warning("6: " + StringUtils::toString(table.cdf(xv, Yv, 100, 1e-8, 1e-4, 0.0004, 1e7), 'e', 9));
        log()->warning("7: " + std::to_string(xv.size()));

        TextOutFile out(this, "stellar_sed", "stellar sed");
        out.addColumn("wavelength");
        out.addColumn("normalized cumulative SED");
        for (size_t i = 0; i!=xv.size(); ++i) out.writeRow(xv[i],Yv[i]);
    }
}

////////////////////////////////////////////////////////////////////
