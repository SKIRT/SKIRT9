/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MonteCarloSimulation.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PhotonPacket.hpp"
#include "ProcessManager.hpp"
#include "SecondarySourceSystem.hpp"
#include "ShortArray.hpp"
#include "SpecialFunctions.hpp"
#include "StringUtils.hpp"
#include "TimeLogger.hpp"

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setupSimulation()
{
    // perform regular setup for the hierarchy and wait for all processes to finish
    {
        TimeLogger logger(log(), "setup");
        _config->setup();  // first of all perform setup for the configuration object
        SimulationItem::setup();
        wait("setup");
    }

    // write setup output
    {
        TimeLogger logger(log(), "setup output");

        // notify the probe system
        probeSystem()->probeSetup();
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setupSelfBefore()
{
    Simulation::setupSelfBefore();

    // construct a secondary source system to help launch secondary photon packets if required
    if (_config->hasSecondaryEmission()) _secondarySourceSystem = new SecondarySourceSystem(this);
}

////////////////////////////////////////////////////////////////////

Configuration* MonteCarloSimulation::config() const
{
    return _config;
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runSimulation()
{
    // run the simulation
    {
        TimeLogger logger(log(), "the run");

        bool hasPrimaryLuminosity = sourceSystem()->luminosity() > 0.;

        // special case of merged primary and secondary iterations
        if (_config->hasMergedIterations() && hasPrimaryLuminosity)
        {
            runPrimaryEmissionIterations();
            runMergedEmissionIterations();
            runPrimaryEmission();
            runSecondaryEmission();
        }
        else
        {
            // primary emission phase, possibly with dynamic medium state iterations
            if (_config->hasPrimaryIterations() && hasPrimaryLuminosity) runPrimaryEmissionIterations();
            runPrimaryEmission();

            // optional secondary emission phase, possibly with dynamic secondary emission iterations
            if (_config->hasSecondaryEmission())
            {
                if (_config->hasSecondaryIterations()) runSecondaryEmissionIterations();
                runSecondaryEmission();
            }
        }
    }

    // write final output
    {
        TimeLogger logger(log(), "final output");

        // notify the probe system
        probeSystem()->probeRun();

        // write instrument output
        instrumentSystem()->flush();
        instrumentSystem()->write();
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runPrimaryEmission()
{
    string segment = "primary emission";
    TimeLogger logger(log(), segment);

    // clear the radiation field
    if (_config->hasRadiationField()) mediumSystem()->clearRadiationField(true);

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
        initProgress(segment, Npp);
        sourceSystem()->prepareForLaunch(Npp);
        auto parallel = find<ParallelFactory>()->parallelDistributed();
        parallel->call(
            Npp, [this](size_t i, size_t n) { performLifeCycle(i, n, true, true, _config->hasRadiationField()); });
        instrumentSystem()->flush();
    }

    // wait for all processes to finish and synchronize the radiation field
    wait(segment);
    if (_config->hasRadiationField()) mediumSystem()->communicateRadiationField(true);

    // update semi-dynamic medium state if needed
    if (_config->hasSemiDynamicState()) mediumSystem()->updateSemiDynamicMediumState();
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runSecondaryEmission()
{
    string segment = "secondary emission";
    TimeLogger logger(log(), segment);

    // determine whether we need to store the radiation field during secondary emission
    // if so, clear the secondary radiation field
    bool storeRF = _config->storeEmissionRadiationField();
    if (storeRF) mediumSystem()->clearRadiationField(false);

    // shoot photons from secondary sources, if needed
    size_t Npp = _config->numSecondaryPackets();
    if (!Npp)
    {
        log()->warning("Skipping secondary emission because no photon packets were requested");
    }
    else if (!_secondarySourceSystem->prepareForLaunch(Npp))
    {
        log()->warning("Skipping secondary emission because the total luminosity of secondary sources is zero");
    }
    else
    {
        initProgress(segment, Npp);
        auto parallel = find<ParallelFactory>()->parallelDistributed();
        parallel->call(Npp, [this, storeRF](size_t i, size_t n) { performLifeCycle(i, n, false, true, storeRF); });
        instrumentSystem()->flush();
    }

    // wait for all processes to finish and synchronize the radiation field if needed
    wait(segment);
    if (storeRF) mediumSystem()->communicateRadiationField(false);
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runPrimaryEmissionIterations()
{
    // when this function is called
    //  - the number of photon packets and the source luminosity are guaranteed to be nonzero
    //  - the data structures to store the radiation field are guaranteed to exist

    // get the parallel engine
    auto parallel = find<ParallelFactory>()->parallelDistributed();

    // get the parameters controlling the dynamic state iteration
    size_t Npp = _config->numPrimaryIterationPackets();
    int minIters = _config->minPrimaryIterations();
    int maxIters = _config->maxPrimaryIterations();

    // prepare the source system for the appropriate number of packets
    sourceSystem()->prepareForLaunch(Npp);

    // loop over the dynamic state iterations; the loop exits
    //   - if convergence is reached after the minimum number of iterations, or
    //   - if the maximum number of iterations has completed, even if there is no convergence
    int iter = 0;
    while (true)
    {
        ++iter;
        bool converged = false;
        {
            string segment = "primary emission iteration " + std::to_string(iter);
            TimeLogger logger(log(), segment);

            // clear the radiation field
            mediumSystem()->clearRadiationField(true);

            // launch photon packets
            initProgress(segment, Npp);
            parallel->call(Npp, [this](size_t i, size_t n) { performLifeCycle(i, n, true, false, true); });
            instrumentSystem()->flush();

            // wait for all processes to finish and synchronize the radiation field
            wait(segment);
            mediumSystem()->communicateRadiationField(true);

            // update the medium state based on the newly established radiation field
            converged = mediumSystem()->updateDynamicMediumState();
        }

        // force at least the minimum number of iterations
        if (converged && iter < minIters)
        {
            log()->info("Convergence reached but continuing until " + std::to_string(minIters)
                        + " iterations have been performed");
        }
        // exit the loop if convergence has been reached after at least the minimum number of iterations
        else if (converged && iter >= minIters)
        {
            log()->info("Convergence reached after " + std::to_string(iter) + " iterations");
            break;
        }
        // continue if convergence has not been reached after fewer than the maximum number of iterations
        else if (!converged && iter < maxIters)
        {
            log()->info("Convergence not yet reached after " + std::to_string(iter) + " iterations");
        }
        // exit the loop if convergence has not been reached after the maximum number of iterations
        else if (!converged && iter >= maxIters)
        {
            log()->error("Convergence not yet reached after " + std::to_string(iter) + " iterations");
            break;
        }
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runSecondaryEmissionIterations()
{
    // get number of photons (guaranteed to be nonzero)
    size_t Npp = _config->numSecondaryIterationPackets();

    // get the parallel engine
    auto parallel = find<ParallelFactory>()->parallelDistributed();

    // get the parameters controlling the secondary iteration
    int minIters = _config->minSecondaryIterations();
    int maxIters = _config->maxSecondaryIterations();
    double fractionOfPrimary = _config->maxFractionOfPrimary();
    double fractionOfPrevious = _config->maxFractionOfPrevious();

    // initialize the total absorbed luminosity in the previous iteration
    double prevLabsdust = 0.;

    // iterate over the maximum number of iterations; the loop body returns from the function
    // when convergence is reached after the minimum number of iterations have been completed
    for (int iter = 1; iter <= maxIters; iter++)
    {
        string segment = "secondary emission iteration " + std::to_string(iter);
        {
            TimeLogger logger(log(), segment);

            // clear the secondary radiation field
            mediumSystem()->clearRadiationField(false);

            // prepare the source system; terminate if secondary luminosity is zero (which would be very unusual)
            if (!_secondarySourceSystem->prepareForLaunch(Npp))
            {
                log()->warning(
                    "Skipping secondary emission iterations because the total luminosity of secondary sources is zero");
                return;
            }

            // launch photon packets
            initProgress(segment, Npp);
            parallel->call(Npp, [this](size_t i, size_t n) { performLifeCycle(i, n, false, false, true); });
            instrumentSystem()->flush();

            // wait for all processes to finish and synchronize the radiation field
            wait(segment);
            mediumSystem()->communicateRadiationField(false);

            // update semi-dynamic medium state if needed
            if (_config->hasSemiDynamicState()) mediumSystem()->updateSemiDynamicMediumState();
        }

        // determine and log the total absorbed luminosity
        double Labsprim, Labsseco;
        std::tie(Labsprim, Labsseco) = mediumSystem()->totalDustAbsorbedLuminosity();
        log()->info("The total dust-absorbed primary luminosity is "
                    + StringUtils::toString(units()->obolluminosity(Labsprim), 'g') + " " + units()->ubolluminosity());
        log()->info("The total dust-absorbed secondary luminosity in iteration " + std::to_string(iter) + " is "
                    + StringUtils::toString(units()->obolluminosity(Labsseco), 'g') + " " + units()->ubolluminosity());

        // log the current performance and corresponding convergence criteria
        if (Labsprim > 0. && Labsseco > 0.)
        {
            if (iter == 1)
            {
                log()->info("--> absorbed secondary luminosity is "
                            + StringUtils::toString(Labsseco / Labsprim * 100., 'f', 2)
                            + "% of absorbed primary luminosity (convergence criterion is "
                            + StringUtils::toString(fractionOfPrimary * 100., 'f', 2) + "%)");
            }
            else
            {
                log()->info("--> absorbed secondary luminosity changed by "
                            + StringUtils::toString(abs((Labsseco - prevLabsdust) / Labsseco) * 100., 'f', 2)
                            + "% compared to previous iteration (convergence criterion is "
                            + StringUtils::toString(fractionOfPrevious * 100., 'f', 2) + "%)");
            }
        }

        // force at least the minimum number of iterations
        if (iter < minIters)
        {
            log()->info("Continuing until " + std::to_string(minIters) + " iterations have been performed");
        }
        else
        {
            // the self-absorption iteration has reached convergence if one or more of the following conditions holds:
            // - the absorbed primary luminosity is zero
            // - the absorbed secondary luminosity is zero
            // - the absorbed secondary luminosity is less than a given fraction of the absorbed primary luminosity
            // - the absorbed secondary luminosity has changed by less than a given fraction compared to previous iter
            if (Labsprim <= 0. || Labsseco <= 0. || Labsseco / Labsprim < fractionOfPrimary
                || abs((Labsseco - prevLabsdust) / Labsseco) < fractionOfPrevious)
            {
                log()->info("Convergence reached after " + std::to_string(iter) + " iterations");
                return;  // end the iteration by returning from the function
            }
            else
            {
                log()->info("Convergence not yet reached after " + std::to_string(iter) + " iterations");
            }
        }
        prevLabsdust = Labsseco;
    }

    // if the loop runs out, convergence was not reached even after the maximum number of iterations
    log()->error("Convergence not yet reached after " + std::to_string(maxIters) + " iterations");
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runMergedEmissionIterations()
{
    throw FATALERROR("Including primary emission in secondary emission iterations is not yet implemented");
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

    log()->info("Launching " + StringUtils::toString(static_cast<double>(numTotal)) + " " + _segment
                + " photon packets");
    log()->infoSetElapsed(numTotal);
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::logProgress(size_t numDone)
{
    // log message if the minimum time has elapsed
    log()->infoIfElapsed("Launched " + _segment + " photon packets: ", numDone);
}

////////////////////////////////////////////////////////////////////

namespace
{
    // maximum number of photon packets processed between two invocations of infoIfElapsed()
    const size_t logProgressChunkSize = 1000;
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::performLifeCycle(size_t firstIndex, size_t numIndices, bool primary, bool peel, bool store)
{
    PhotonPacket pp, ppp;

    // loop over the history indices, with interruptions for progress logging
    while (numIndices)
    {
        size_t currentChunkSize = min(logProgressChunkSize, numIndices);
        for (size_t historyIndex = firstIndex; historyIndex != firstIndex + currentChunkSize; ++historyIndex)
        {
            // launch a photon packet from the requested source
            if (primary)
                sourceSystem()->launch(&pp, historyIndex);
            else
                _secondarySourceSystem->launch(&pp, historyIndex);
            if (pp.luminosity() > 0)
            {
                if (peel) peelOffEmission(&pp, &ppp);

                // trace the packet through the media, if any
                if (_config->hasMedium())
                {
                    if (_config->forceScattering())
                    {
                        // perform cycle with forced scattering
                        double Lthreshold = pp.luminosity() / _config->minWeightReduction();
                        int minScattEvents = _config->minScattEvents();
                        while (true)
                        {
                            // calculate segments and optical depths for the complete path
                            mediumSystem()->setOpticalDepths(&pp);

                            // advance the packet
                            if (store) storeRadiationField(&pp);
                            simulateForcedPropagation(&pp);

                            // if the packet's weight drops below the threshold, terminate it
                            if (pp.luminosity() <= 0
                                || (pp.luminosity() <= Lthreshold && pp.numScatt() >= minScattEvents))
                                break;

                            // process the scattering event
                            if (peel) peelOffScattering(&pp, &ppp);
                            mediumSystem()->simulateScattering(random(), &pp);
                        }
                    }
                    else
                    {
                        // perform cycle without forced scattering
                        while (true)
                        {
                            // generate a random interaction optical depth
                            double tauscat = random()->expon();

                            // find the physical interaction point corresponding to this optical depth
                            // if the interaction optical depth is outside of the path, terminate the photon packet
                            if (!mediumSystem()->setInteractionPoint(&pp, tauscat)) break;

                            // advance the packet
                            simulateNonForcedPropagation(&pp);

                            // process the scattering event
                            if (peel) peelOffScattering(&pp, &ppp);
                            mediumSystem()->simulateScattering(random(), &pp);
                        }
                    }
                }
            }
        }

        // log progress
        logProgress(currentChunkSize);
        firstIndex += currentChunkSize;
        numIndices -= currentChunkSize;
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::peelOffEmission(const PhotonPacket* pp, PhotonPacket* ppp)
{
    for (Instrument* instrument : _instrumentSystem->instruments())
    {
        if (!instrument->isSameObserverAsPreceding())
        {
            const Direction bfkobs = instrument->bfkobs(pp->position());
            ppp->launchEmissionPeelOff(pp, bfkobs);

            // if the photon packet is polarised, we have to rotate the Stokes vector into the frame of the instrument
            if (ppp->isPolarized())
            {
                ppp->rotateIntoPlane(bfkobs, instrument->bfky(pp->position()));
            }
        }
        instrument->detect(ppp);
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::storeRadiationField(const PhotonPacket* pp)
{
    // use a faster version in case there are no kinematics
    if (_config->hasConstantPerceivedWavelength())
    {
        int ell = _config->radiationFieldWLG()->bin(pp->wavelength());
        if (ell >= 0)
        {
            double luminosity = pp->luminosity();
            bool hasPrimaryOrigin = pp->hasPrimaryOrigin();

            double lnExtBeg = 0.;  // extinction factor and its logarithm at begin of current segment
            double extBeg = 1.;
            for (const auto& segment : pp->segments())
            {
                double lnExtEnd = -segment.tau;  // extinction factor and its logarithm at end of current segment
                double extEnd = exp(lnExtEnd);
                int m = segment.m;
                if (m >= 0)
                {
                    // use this flavor of the lnmean function to avoid recalculating the logarithm of the extinction
                    double extMean = SpecialFunctions::lnmean(extEnd, extBeg, lnExtEnd, lnExtBeg);
                    double Lds = luminosity * extMean * segment.ds;
                    mediumSystem()->storeRadiationField(hasPrimaryOrigin, m, ell, Lds);
                }
                lnExtBeg = lnExtEnd;
                extBeg = extEnd;
            }
        }
    }
    else
    {
        double lnExtBeg = 0.;  // extinction factor and its logarithm at begin of current segment
        double extBeg = 1.;
        for (const auto& segment : pp->segments())
        {
            double lnExtEnd = -segment.tau;  // extinction factor and its logarithm at end of current segment
            double extEnd = exp(lnExtEnd);
            int m = segment.m;
            if (m >= 0)
            {
                double lambda = pp->perceivedWavelength(mediumSystem()->bulkVelocity(m),
                                                        _config->hubbleExpansionRate() * segment.s);
                int ell = _config->radiationFieldWLG()->bin(lambda);
                if (ell >= 0)
                {
                    // use this flavor of the lnmean function to avoid recalculating the logarithm of the extinction
                    double extMean = SpecialFunctions::lnmean(extEnd, extBeg, lnExtEnd, lnExtBeg);
                    double Lds = pp->perceivedLuminosity(lambda) * extMean * segment.ds;
                    mediumSystem()->storeRadiationField(pp->hasPrimaryOrigin(), m, ell, Lds);
                }
            }
            lnExtBeg = lnExtEnd;
            extBeg = extEnd;
        }
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::simulateForcedPropagation(PhotonPacket* pp)
{
    // get the total optical depth
    double taupath = pp->totalOpticalDepth();

    // if there is no extinction along the path, this photon packet cannot scatter, so terminate it right away
    if (taupath <= 0.)
    {
        pp->applyBias(0.);
        return;
    }

    // generate a random optical depth
    double xi = _config->pathLengthBias();
    double tau = 0.;
    if (xi == 0.)
    {
        tau = random()->exponCutoff(taupath);
    }
    else
    {
        tau = random()->uniform() < xi ? random()->uniform() * taupath : random()->exponCutoff(taupath);
        double p = -exp(-tau) / expm1(-taupath);
        double q = (1.0 - xi) * p + xi / taupath;
        double weight = p / q;
        pp->applyBias(weight);
    }

    // determine the physical position of the interaction point
    pp->findInteractionPoint(tau);

    // calculate the albedo for the cell containing the interaction point
    double albedo = mediumSystem()->albedoForScattering(pp);

    // adjust the weight by the scattered fraction
    pp->applyBias(-expm1(-taupath) * albedo);

    // advance the position
    pp->propagate(pp->interactionDistance());
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::simulateNonForcedPropagation(PhotonPacket* pp)
{
    // calculate the albedo for the cell containing the interaction point
    double albedo = mediumSystem()->albedoForScattering(pp);

    // adjust the weight by the albedo
    pp->applyBias(albedo);

    // advance the photon packet position
    pp->propagate(pp->interactionDistance());
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::peelOffScattering(PhotonPacket* pp, PhotonPacket* ppp)
{
    // determine the perceived wavelength at the scattering location
    double lambda = mediumSystem()->perceivedWavelengthForScattering(pp);

    // determine the scattering opacity weight for each medium component;
    // abort if none of the media scatter this photon packet
    ShortArray wv;
    if (!mediumSystem()->weightsForScattering(wv, lambda, pp)) return;

    // now do the actual peel-off for each instrument
    for (Instrument* instr : _instrumentSystem->instruments())
    {
        if (!instr->isSameObserverAsPreceding())
        {
            // get the direction towards the instrument and (for polarization only) its Y-axis orientation
            Direction bfkobs = instr->bfkobs(pp->position());
            Direction bfky = _config->hasPolarization() ? instr->bfky(pp->position()) : Direction();

            // calculate peel-off for all medium components and launch the peel-off photon packet
            // (all media must either support polarization or not; combining these support levels is not allowed)
            mediumSystem()->peelOffScattering(lambda, wv, bfkobs, bfky, pp, ppp);
        }

        // have the peel-off photon packet detected
        instr->detect(ppp);
    }
}

////////////////////////////////////////////////////////////////////
