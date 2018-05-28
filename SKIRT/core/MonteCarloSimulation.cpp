/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MonteCarloSimulation.hpp"
#include "Constants.hpp"
#include "FITSInOut.hpp"
#include "FatalError.hpp"
#include "LockFree.hpp"
#include "Log.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PhotonPacket.hpp"
#include "ProcessManager.hpp"
#include "Random.hpp"
#include "StopWatch.hpp"
#include "StoredTable.hpp"
#include "StringUtils.hpp"
#include "Table.hpp"
#include "TextOutFile.hpp"
#include "TimeLogger.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"
#include <array>
#include <atomic>
#include <map>
#include <mutex>
#include <thread>

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setupSimulation()
{
    // perform regular setup for the hierarchy and wait for all processes to finish
    TimeLogger logger(log(), "setup");
    SimulationItem::setup();
    wait("setup");

    // notify the probe system
    probeSystem()->probeSetup();
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setupSelfBefore()
{
    Simulation::setupSelfBefore();
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setEmulationMode()
{
    _emulationMode = true;
    _numPackets = 0;
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

void MonteCarloSimulation::wait(std::string scope)
{
    if (ProcessManager::isMultiProc())
    {
        find<Log>()->info("Waiting for other processes to finish " + scope + "...");
        ProcessManager::wait();
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runSimulation()
{
    {
        // perform the run and wait for all processes to finish
        TimeLogger logger(log(), "the run");
        test();
        wait("the run");

        // notify the probe system
        probeSystem()->probeRun();
    }
    {
        // write simulation output
        TimeLogger logger(log(), "the output");
        instrumentSystem()->flush();
        instrumentSystem()->write();
    }
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

void MonteCarloSimulation::doSolarEmissionChunk(size_t firstIndex, size_t numIndices)
{
    PhotonPacket pp,ppp;
    double Lpacket = _Ltot / numPackets();
    double sigma = 0.01 * Constants::AU();

    // iterate over the chunk of indices
    for (size_t historyIndex = firstIndex; historyIndex!=firstIndex+numIndices; ++historyIndex)
    {
        // generate a random position for this source according to a Gaussian kernel
        double x = random()->gauss();
        double y = random()->gauss();
        double z = random()->gauss();
        Position bfr( Vec(x,y,z)*sigma );

        // generate a random wavelength for this source from the SED
        double lambda = random()->cdf(_sunLambda, _sunCDF);

        // launch the primary photon packet
        pp.launch(historyIndex, lambda, Lpacket, bfr, random()->direction());
        pp.setPrimaryOrigin(0);

        // peel off towards each instrument
        for (Instrument* instrument : _instrumentSystem->instruments())
        {
            ppp.launchEmissionPeelOff(&pp, instrument->bfkobs(pp.position()), 1.);
            instrument->detect(&ppp);
        }
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::test()
{
    // --- parallelization ---
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

    // --- sampling SEDs and calibrating instruments ---
    {
        // get the cumulative distribution for the solar SED
        StoredTable<1> table;
        table.open(this, "SunSED", "lambda(m)", "Llambda(W/m)");
        _Ltot = table.cdf(_sunLambda, _sunCDF, 200, 0.1e-6, 20e-6);
        log()->info("Fraction of solar luminosity: " + StringUtils::toString(_Ltot/Constants::Lsun() , 'f', 5));

        // shoot photons
        auto parallel = find<ParallelFactory>()->parallelDistributed();
        parallel->call([this](size_t i ,size_t n) { doSolarEmissionChunk(i, n); }, numPackets());
        instrumentSystem()->flush();
    }
}

////////////////////////////////////////////////////////////////////
