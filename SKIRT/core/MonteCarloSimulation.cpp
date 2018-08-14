/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MonteCarloSimulation.hpp"
#include "Log.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PhotonPacket.hpp"
#include "ProcessManager.hpp"
#include "StopWatch.hpp"
#include "StringUtils.hpp"
#include "TimeLogger.hpp"

// uncomment this line to produce chunk debugging messages
//#define LOG_CHUNKS

#ifdef LOG_CHUNKS
#include <map>
#include <mutex>
#include <thread>
#endif

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setupSimulation()
{
    // perform regular setup for the hierarchy and wait for all processes to finish
    TimeLogger logger(log(), "setup");
    _config->setup();           // first of all perform setup for the configuration object
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

void MonteCarloSimulation::setupSelfAfter()
{
    Simulation::setupSelfAfter();

    log()->info("Model symmetry is " + std::to_string(dimension()) + "D");
}

////////////////////////////////////////////////////////////////////

Configuration* MonteCarloSimulation::config() const
{
    return _config;
}

////////////////////////////////////////////////////////////////////

int MonteCarloSimulation::dimension() const
{
    return sourceSystem()->dimension();
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runSimulation()
{
    {
        TimeLogger logger(log(), "the run");

        // shoot photons from primary sources, if needed
        size_t Npp = _config->numPrimaryPackets();
        if (!Npp)
        {
            log()->warning("Skipping primary emission because no photon packets were requested");
        }
        else if (!sourceSystem()->luminosity())
        {
            log()->warning("Skipping primary emission because the total luminosity of primary sources is zero");
        }
        else
        {
            initProgress("primary emission", Npp);
            sourceSystem()->prepareForlaunch(Npp);
            auto parallel = find<ParallelFactory>()->parallelDistributed();
            StopWatch::start();
            parallel->call(Npp, [this](size_t i ,size_t n) { doPrimaryEmissionChunk(i, n); });
            StopWatch::stop();
            instrumentSystem()->flush();
        }

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

void MonteCarloSimulation::wait(std::string scope)
{
    if (ProcessManager::isMultiProc())
    {
        log()->info("Waiting for other processes to finish " + scope + "...");
        ProcessManager::wait();
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::initProgress(string segment, size_t numTotal)
{
    _segment = segment;
    _numTotal = numTotal;

    log()->info("Launching " + StringUtils::toString(static_cast<double>(_numTotal))
                + " " + _segment + " photon packages");
    log()->infoSetElapsed();
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::logProgress(size_t numDone)
{
    // log message if the minimum time has elapsed
    log()->infoIfElapsed("Launched " + _segment + " photon packages: ", numDone, _numTotal);
}

////////////////////////////////////////////////////////////////////

namespace
{
    // maximum number of photon packets processed between two invocations of logProgress()
    const size_t logProgressChunkSize = 50000;
}

////////////////////////////////////////////////////////////////////

#ifdef LOG_CHUNKS
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
#endif

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::doPrimaryEmissionChunk(size_t firstIndex, size_t numIndices)
{
#ifdef LOG_CHUNKS
    // log chunk info for debugging purposes
    int id = threadID();
    log()->warning("[T" + std::to_string(id) + "] Chunk: "
                   + std::to_string(firstIndex) + "," + std::to_string(numIndices));

    // time one of the threads for debugging purposes
    if (id==1) StopWatch::start();
#endif

    // actually shoot the photon packets
    {
        PhotonPacket pp,ppp;

        size_t endIndex = firstIndex + numIndices;
        for (size_t historyIndex=firstIndex; historyIndex!=endIndex; ++historyIndex)
        {
            if (historyIndex%logProgressChunkSize == 0) logProgress(historyIndex);

            // launch a photon packet
            sourceSystem()->launch(&pp, historyIndex);

            // peel off towards each instrument
            for (Instrument* instrument : _instrumentSystem->instruments())
            {
                ppp.launchEmissionPeelOff(&pp, instrument->bfkobs(pp.position()));
                instrument->detect(&ppp);
            }
        }
        logProgress(endIndex);
    }

#ifdef LOG_CHUNKS
    // time one of the threads for debugging purposes
    if (id==1) StopWatch::stop();
#endif
}

////////////////////////////////////////////////////////////////////
