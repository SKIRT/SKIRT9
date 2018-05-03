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
#include "ProcessManager.hpp"
#include "StopWatch.hpp"
#include "StoredTable.hpp"
#include "Table.hpp"
#include "TextOutFile.hpp"
#include <atomic>
#include <map>
#include <mutex>
#include <thread>

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
    // get a short but consistent identifier for the current thread
    int threadID()
    {
        static std::map<std::thread::id, int> threads;
        static std::mutex mutex;
        auto me = std::this_thread::get_id();
        std::unique_lock<std::mutex> lock(mutex);
        if (!threads.count(me)) threads.emplace(me, threads.size()+1);
        return threads.at(me);
    }

    // test function adds 1/inverseFraction to numPixels random pixels in the specified frame
    void addPixels(Table<2>& frame, Random* random, size_t inverseFraction, size_t firstIndex, size_t numIndices)
    {
        int id = threadID();
        random->find<Log>()->warning("[T" + std::to_string(id) + "] Chunk: "
                                     + std::to_string(firstIndex) + "," + std::to_string(numIndices));

        ///if (firstIndex>10*1000*1000) throw FATALERROR("Test exception handling");

        if (id==1) StopWatch::start();      // time only one of the threads
        for (size_t p = 0; p!=numIndices; ++p)
        {
            size_t i = static_cast<size_t>( random->uniform() * frame.size(0) );
            size_t j = static_cast<size_t>( random->uniform() * frame.size(1) );
            LockFree::add(frame(i,j), 1./inverseFraction);
        }
        if (id==1) StopWatch::stop();
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runSelf()
{
    TimeLogger logger(log(), "the test phase");

    {
        StoredTable<1> table;
        table.open(this, "SunSED", "lambda(m)", "Llambda(W/m)");

        log()->info("0: " + StringUtils::toString(table[0.55e-6], 'e', 9));
        log()->info("1: " + StringUtils::toString(table[1e-6], 'e', 9));
        log()->info("2: " + StringUtils::toString(table[100e-6], 'e', 9));
    }
    {
        StoredTable<2> table;
        table.open(this, "DustEM_Gra_OpticalProps", "lambda(m),a(m)", "Qabs(1)");

        log()->info("3: " + StringUtils::toString(table(0.55e-6, 0.01e-6), 'e', 9));
    }
    {
        StoredTable<1> table;
        table.open(this, "MeanIvezicBenchmarkOpticalProps", "lambda(m)", "sigmaabs(m2/H)");

        Array xv, Yv;
        log()->info("4: " + StringUtils::toString(table.cdf(xv, Yv, 100, 0.1e-6, 1e-6), 'e', 9));
        log()->info("5: " + StringUtils::toString(table.cdf(xv, Yv, 1000, 1e-6, 10e-6), 'e', 9));
    }
    {
        StoredTable<3> table;
        table.open(this, "BruzualCharlotSEDFamily_Chabrier_hr", "lambda(m),Z(1),t(yr)", "Llambda(W/m)");

        Array xv, Yv;
        log()->info("6: " + StringUtils::toString(table.cdf(xv, Yv, 100, 1e-8, 1e-4, 0.0004, 1e7), 'e', 9));
        log()->info("7: " + std::to_string(xv.size()));

        TextOutFile out(this, "stellar_sed", "stellar sed");
        out.addColumn("wavelength");
        out.addColumn("normalized cumulative SED");
        for (size_t i = 0; i!=xv.size(); ++i) out.writeRow(xv[i],Yv[i]);
    }
    {
        const size_t numPixels = 100*1000*1000 + 11;
        Table<2> frame(500,500);

        auto parallel = find<ParallelFactory>()->parallelDistributed();
        StopWatch::start();
        parallel->call([this,&frame](size_t i ,size_t n) { addPixels(frame, random(), numPixels, i, n); }, numPixels);
        StopWatch::stop();
        ProcessManager::sumToRoot(frame.data());

        if (ProcessManager::isRoot())
            log()->info("Frame intensity: " + StringUtils::toString(frame.data().sum(), 'e', 9));

        FITSInOut::write(this, "frame", "frame", frame.data(), frame.size(0), frame.size(1), 1, 0,0,0,0,"","");
    }
}

////////////////////////////////////////////////////////////////////
