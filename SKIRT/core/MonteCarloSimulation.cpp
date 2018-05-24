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

/*
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
*/

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runSelf()
{
    TimeLogger logger(log(), "the test phase");
/*
    // --- stored tables ---
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

    // --- parallelization ---
    {
        const size_t numPixels = 1*1000*1000 + 11;
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

    // --- wavelength grids ---
    {
        WavelengthGrid* lambdagrid = instrumentSystem()->defaultWavelengthGrid();
        if (lambdagrid)
        {
            Units* units = find<Units>();

            // create a text file and add the columns
            TextOutFile file(this, "wavelengths", "wavelengths");
            file.addColumn("characteristic wavelength (" + units->uwavelength() + ")", 'e', 9);
            file.addColumn("wavelength bin width (" + units->uwavelength() + ")", 'e', 9);
            file.addColumn("left border of wavelength bin (" + units->uwavelength() + ")", 'e', 9);
            file.addColumn("right border of wavelength bin (" + units->uwavelength() + ")", 'e', 9);

            // write the rows
            int Nlambda = lambdagrid->numWavelengths();
            for (int ell=0; ell!=Nlambda; ++ell)
            {
                file.writeRow(units->owavelength(lambdagrid->lambda(ell)),
                              units->owavelength(lambdagrid->dlambda(ell)),
                              units->owavelength(lambdagrid->lambdaLeft(ell)),
                              units->owavelength(lambdagrid->lambdaRight(ell)) );
            }
        }
    }
    {
        WavelengthGrid* lambdagrid = instrumentSystem()->defaultWavelengthGrid();
        if (lambdagrid)
        {
            std::vector<double> wavelengths({0.5, 1, 1.5, 2, 2.05, 7});  // in micron
            for (double wavelength : wavelengths)
            {
                log()->info("Wavelength: " + StringUtils::toString(wavelength, 'f', 2)
                            + "  Index: " + std::to_string(lambdagrid->ell(wavelength*1e-6)));
            }
        }
    }

    // --- instruments ---
    {
        auto parallel = find<ParallelFactory>()->parallelDistributed();
        parallel->call([this](size_t i ,size_t n) { doTestEmissionChunk(i, n); }, numPackets());
        instrumentSystem()->flush();
        instrumentSystem()->write();
    }
*/
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
        instrumentSystem()->write();
    }
}

////////////////////////////////////////////////////////////////////

/*
namespace
{
    const int numSources = 2;
    const std::array<double,numSources> sourceSize{{30, 80}};                    //  in parsec
    const std::array<Vec,numSources> sourceCenter{{Vec(100,0,0), Vec(0,200,0)}}; //  in parsec
    const std::array<double,numSources> sourceLuminosity{{1, 10}};               //  in Lsun
    const std::array<double,numSources> sourceMinLambda{{0.15, 4}};              //  in micron
    const std::array<double,numSources> sourceMaxLambda{{20, 50}};               //  in micron
}
*/

////////////////////////////////////////////////////////////////////

/*
void MonteCarloSimulation::doTestEmissionChunk(size_t firstIndex, size_t numIndices)
{
    PhotonPacket pp,ppp;
    double micron = 1e-6;
    double pc = Constants::pc();
    double LsunPerPacket = Constants::Lsun() / numPackets();

    // iterate over the chunk of indices
    for (size_t historyIndex = firstIndex; historyIndex!=firstIndex+numIndices; ++historyIndex)
    {
        // randomly select one of the sources
        int index = random()->uniform() * numSources;

        // generate a random position for this source according to a Gaussian kernel
        double x = random()->gauss();
        double y = random()->gauss();
        double z = random()->gauss();
        Position bfr( (sourceCenter[index] + Vec(x,y,z)*sourceSize[index]) * pc );

        // generate a random wavelength for this source, uniform in the given wavelength interval
        double lambda = sourceMinLambda[index] + random()->uniform() * (sourceMaxLambda[index]-sourceMinLambda[index]);
        lambda *= micron;

        // get the equal-share luminosity
        double luminosity = sourceLuminosity[index] * LsunPerPacket;

        // launch the primary photon packet
        pp.launch(historyIndex, lambda, luminosity, bfr, random()->direction());
        pp.setPrimaryOrigin(0);

        // peel off towards each instrument
        for (Instrument* instrument : _instrumentSystem->instruments())
        {
            ppp.launchEmissionPeelOff(&pp, instrument->bfkobs(pp.position()), 1.);
            instrument->detect(&ppp);
        }
    }
}
*/

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
