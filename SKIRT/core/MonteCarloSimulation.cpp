/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MonteCarloSimulation.hpp"
#include "Constants.hpp"
#include "FITSInOut.hpp"
#include "FatalError.hpp"
#include "GeometricSource.hpp"
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
        TimeLogger logger(log(), "the run");

        // shoot photons from primary sources
        size_t Npp = numPackets() * sourceSystem()->emissionMultiplier();
        sourceSystem()->prepareForlaunch(Npp);
        auto parallel = find<ParallelFactory>()->parallelDistributed();
        StopWatch::start();
        parallel->call([this](size_t i ,size_t n) { doPrimaryEmissionChunk(i, n); }, Npp);
        StopWatch::stop();
        instrumentSystem()->flush();

        // wait for all processes to finish
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
    // get a short but consistent identifier for the current thread, used for debugging purposes
    int threadID()
    {
        static std::map<std::thread::id, int> threads;
        static std::mutex mutex;
        auto me = std::this_thread::get_id();
        std::unique_lock<std::mutex> lock(mutex);
        if (!threads.count(me)) threads.emplace(me, threads.size()+1);
        return threads.at(me);
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::doPrimaryEmissionChunk(size_t firstIndex, size_t numIndices)
{
    // log chunk info for debugging purposes
    int id = threadID();
    log()->warning("[T" + std::to_string(id) + "] Chunk: "
                   + std::to_string(firstIndex) + "," + std::to_string(numIndices));

    // time one of the threads for debugging purposes
    if (id==1) StopWatch::start();

    // actually shoot the photon packets
    {
        PhotonPacket pp,ppp;

        // iterate over the chunk of indices
        for (size_t historyIndex = firstIndex; historyIndex!=firstIndex+numIndices; ++historyIndex)
        {
            // launch a photon packet
            sourceSystem()->launch(&pp, historyIndex);

            // peel off towards each instrument
            for (Instrument* instrument : _instrumentSystem->instruments())
            {
                ppp.launchEmissionPeelOff(&pp, instrument->bfkobs(pp.position()), 1.);
                instrument->detect(&ppp);
            }
        }
    }

    // time one of the threads for debugging purposes
    if (id==1) StopWatch::stop();
}

////////////////////////////////////////////////////////////////////
