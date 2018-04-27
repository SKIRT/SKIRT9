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
#include "FITSInOut.hpp"
#include "LockFree.hpp"
#include "ParallelFactory.hpp"
#include "Parallel.hpp"
#include "StoredTable.hpp"
#include "Table.hpp"
#include "TextOutFile.hpp"
#include <atomic>

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

namespace
{
    // test function adds 1/inverseFraction to numPixels random pixels in the specified frame
    void addPixels(Table<2>& frame, Random* random, size_t inverseFraction, size_t firstIndex, size_t numIndices)
    {
        for (size_t p = 0; p!=numIndices; ++p)
        {
            size_t i = static_cast<size_t>( random->uniform() * frame.size(0) );
            size_t j = static_cast<size_t>( random->uniform() * frame.size(1) );
            LockFree::add(frame(i,j), 1./inverseFraction);
        }
        (void)firstIndex;
    }
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
    {
        const size_t numPixels = 100*1000*1000 + 11;
        Table<2> frame(500,500);

        auto parallel = find<ParallelFactory>()->parallelDistributed();
        parallel->call([this,&frame](size_t i ,size_t n) { addPixels(frame, random(), numPixels, i, n); }, numPixels);
        log()->warning("Frame intensity: " + StringUtils::toString(frame.data().sum(), 'e', 9));

        FITSInOut::write(this, "frame", "frame", frame.data(), frame.size(0), frame.size(1), 1, 0,0,0,0,"","");
    }
}

////////////////////////////////////////////////////////////////////
