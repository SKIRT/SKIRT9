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

void MonteCarloSimulation::setupSelfAfter()
{
    Simulation::setupSelfAfter();

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
            if (_config->hasPrimaryIterations()) runPrimaryEmissionIterations();
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

    // update secondary dynamic medium state if applicable (in which case we have a medium system)
    if (_config->hasSecondaryDynamicState()) mediumSystem()->updateSecondaryDynamicMediumState();
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

namespace
{
    // this helper class checks secondary emission convergence based on the dust-absorbed luminosity
    class DustAbsorptionConvergence
    {
        double _prevLabsseco{0.};  // remembers the absorbed luminosity in the previous iteration

    public:
        // this function determines and logs the total absorbed luminosity and related percentages
        // it returns true if secondary emission can be considered to be converged, false otherwise
        bool logConvergenceInfo(Log* log, Units* units, MediumSystem* mediumSystem, int iter, double fractionOfPrimary,
                                double fractionOfPrevious)
        {
            // determine and log the total absorbed luminosity
            double Labsprim, Labsseco;
            std::tie(Labsprim, Labsseco) = mediumSystem->totalDustAbsorbedLuminosity();
            log->info("The total dust-absorbed primary luminosity is "
                      + StringUtils::toString(units->obolluminosity(Labsprim), 'g') + " " + units->ubolluminosity());
            log->info("The total dust-absorbed secondary luminosity in iteration " + std::to_string(iter) + " is "
                      + StringUtils::toString(units->obolluminosity(Labsseco), 'g') + " " + units->ubolluminosity());

            // log the current performance and corresponding convergence criteria
            if (Labsprim > 0. && Labsseco > 0.)
            {
                if (iter == 1)
                {
                    log->info("--> absorbed secondary luminosity is "
                              + StringUtils::toString(Labsseco / Labsprim * 100., 'f', 2)
                              + "% of absorbed primary luminosity (convergence criterion is "
                              + StringUtils::toString(fractionOfPrimary * 100., 'f', 2) + "%)");
                }
                else
                {
                    log->info("--> absorbed secondary luminosity changed by "
                              + StringUtils::toString(abs((Labsseco - _prevLabsseco) / Labsseco) * 100., 'f', 2)
                              + "% compared to previous iteration (convergence criterion is "
                              + StringUtils::toString(fractionOfPrevious * 100., 'f', 2) + "%)");
                }
            }

            // secondary emission has reached convergence if one or more of the following conditions holds:
            // - the absorbed primary luminosity is zero
            // - the absorbed secondary luminosity is zero
            // - the absorbed secondary luminosity is less than a given fraction of the absorbed primary luminosity
            // - the absorbed secondary luminosity has changed by less than a given fraction compared to previous iter
            bool converged = Labsprim <= 0. || Labsseco <= 0. || Labsseco / Labsprim < fractionOfPrimary
                             || abs((Labsseco - _prevLabsseco) / Labsseco) < fractionOfPrevious;
            _prevLabsseco = Labsseco;
            return converged;
        }
    };

    // this function logs the convergence status and returns true if the loop should exit, false if it should continue
    // specifically, the loop exits
    //   - if convergence is reached after the minimum number of iterations, or
    //   - if the maximum number of iterations has completed, even if there is no convergence
    bool logLoopConvergence(Log* log, bool converged, int iter, int minIters, int maxIters)
    {
        // force at least the minimum number of iterations
        if (converged && iter < minIters)
        {
            log->info("Convergence reached but continuing until " + std::to_string(minIters)
                      + " iterations have been performed");
            return false;
        }
        // exit the loop if convergence has been reached after at least the minimum number of iterations
        if (converged && iter >= minIters)
        {
            log->info("Convergence reached after " + std::to_string(iter) + " iterations");
            return true;
        }
        // continue if convergence has not been reached after fewer than the maximum number of iterations
        if (!converged && iter < maxIters)
        {
            log->info("Convergence not yet reached after " + std::to_string(iter) + " iterations");
            return false;
        }
        // exit the loop if convergence has not been reached after the maximum number of iterations
        if (!converged && iter >= maxIters)
        {
            log->error("Convergence not yet reached after " + std::to_string(iter) + " iterations");
            return true;
        }
        return true;  // the logic can never get here
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runPrimaryEmissionIterations()
{
    // when this function is called
    //  - the number of photon packets and the primary source luminosity are guaranteed to be nonzero
    //  - the data structures to store the radiation field are guaranteed to exist

    // get the parallel engine
    auto parallel = find<ParallelFactory>()->parallelDistributed();

    // get the parameters controlling the dynamic state iteration
    size_t Npp = _config->numPrimaryIterationPackets();
    int minIters = _config->minPrimaryIterations();
    int maxIters = _config->maxPrimaryIterations();

    // prepare the source system for the appropriate number of packets
    sourceSystem()->prepareForLaunch(Npp);

    // loop over the dynamic state iterations
    int iter = 0;
    while (true)
    {
        ++iter;

        bool converged = true;
        {
            string segment = "primary emission iteration " + std::to_string(iter);
            TimeLogger logger(log(), segment);

            mediumSystem()->beginDynamicMediumStateIteration();

            // clear the radiation field
            mediumSystem()->clearRadiationField(true);

            // launch photon packets
            initProgress(segment, Npp);
            parallel->call(Npp, [this](size_t i, size_t n) { performLifeCycle(i, n, true, false, true); });
            instrumentSystem()->flush();

            // wait for all processes to finish and synchronize the radiation field
            wait(segment);
            mediumSystem()->communicateRadiationField(true);

            // update the primary dynamic medium state and log convergence info
            converged = mediumSystem()->updatePrimaryDynamicMediumState();
        }

        // notify the probe system
        probeSystem()->probePrimary(iter);

        // verify and log loop convergence
        if (logLoopConvergence(log(), converged, iter, minIters, maxIters)) break;
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runSecondaryEmissionIterations()
{
    // when this function is called
    //  - the number of photon packets is guaranteed to be nonzero
    //  - the total luminosity of secondary sources may still be zero
    //  - the data structures to store the radiation field are guaranteed to exist

    // get the parallel engine
    auto parallel = find<ParallelFactory>()->parallelDistributed();

    // get the parameters controlling the dynamic secondary emission iteration
    size_t Npp = _config->numSecondaryIterationPackets();
    int minIters = _config->minSecondaryIterations();
    int maxIters = _config->maxSecondaryIterations();
    double fractionOfPrimary = _config->maxFractionOfPrimary();
    double fractionOfPrevious = _config->maxFractionOfPrevious();

    // helper object to verify convergence of secondary emission
    DustAbsorptionConvergence dustConvergence;

    // loop over the secondary emission iterations
    int iter = 0;
    while (true)
    {
        ++iter;

        bool converged = true;
        {
            string segment = "secondary emission iteration " + std::to_string(iter);
            TimeLogger logger(log(), segment);

            mediumSystem()->beginDynamicMediumStateIteration();

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

            // update secondary dynamic medium state and log convergence info
            converged &= mediumSystem()->updateSecondaryDynamicMediumState();

            // log dust emission convergence info
            if (mediumSystem()->hasDust())
                converged &= dustConvergence.logConvergenceInfo(log(), units(), mediumSystem(), iter, fractionOfPrimary,
                                                                fractionOfPrevious);
        }

        // notify the probe system
        probeSystem()->probeSecondary(iter);

        // verify and log loop convergence
        if (logLoopConvergence(log(), converged, iter, minIters, maxIters)) break;
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runMergedEmissionIterations()
{
    // when this function is called
    //  - the number of photon packets and the primary source luminosity are guaranteed to be nonzero
    //  - the total luminosity of secondary sources may still be zero
    //  - the data structures to store the radiation field are guaranteed to exist

    // get the parallel engine
    auto parallel = find<ParallelFactory>()->parallelDistributed();

    // get the parameters controlling the merged iteration
    size_t Npp1 = _config->numPrimaryIterationPackets();
    size_t Npp2 = _config->numSecondaryIterationPackets();
    int minIters = _config->minSecondaryIterations();
    int maxIters = _config->maxSecondaryIterations();
    double fractionOfPrimary = _config->maxFractionOfPrimary();
    double fractionOfPrevious = _config->maxFractionOfPrevious();

    // prepare the primary source system for the appropriate number of packets
    sourceSystem()->prepareForLaunch(Npp1);

    // helper object to verify convergence of secondary emission
    DustAbsorptionConvergence dustConvergence;

    // loop over the merged iterations
    int iter = 0;
    while (true)
    {
        ++iter;

        bool converged = true;
        {
            string segment = "merged primary and secondary emission iteration " + std::to_string(iter);
            string segment1 = "merged primary emission iteration " + std::to_string(iter);
            string segment2 = "merged secondary emission iteration " + std::to_string(iter);
            TimeLogger logger(log(), segment);

            mediumSystem()->beginDynamicMediumStateIteration();

            // clear the radiation field
            mediumSystem()->clearRadiationField(true);

            // launch photon packets
            initProgress(segment1, Npp1);
            parallel->call(Npp1, [this](size_t i, size_t n) { performLifeCycle(i, n, true, false, true); });
            instrumentSystem()->flush();

            // wait for all processes to finish and synchronize the radiation field
            wait(segment1);
            mediumSystem()->communicateRadiationField(true);

            // update secondary dynamic medium state and log convergence info
            converged &= mediumSystem()->updateSecondaryDynamicMediumState();

            // clear the secondary radiation field
            mediumSystem()->clearRadiationField(false);

            // prepare the source system; terminate if secondary luminosity is zero (which would be very unusual)
            if (!_secondarySourceSystem->prepareForLaunch(Npp2))
            {
                log()->warning(
                    "Skipping merged emission iterations because the total luminosity of secondary sources is zero");
                return;
            }

            // launch photon packets
            initProgress(segment2, Npp2);
            parallel->call(Npp2, [this](size_t i, size_t n) { performLifeCycle(i, n, false, false, true); });
            instrumentSystem()->flush();

            // wait for all processes to finish and synchronize the radiation field
            wait(segment2);
            mediumSystem()->communicateRadiationField(false);

            // update the primary dynamic medium state and log convergence info
            converged &= mediumSystem()->updatePrimaryDynamicMediumState();

            // log dust emission convergence info
            if (mediumSystem()->hasDust())
                converged &= dustConvergence.logConvergenceInfo(log(), units(), mediumSystem(), iter, fractionOfPrimary,
                                                                fractionOfPrevious);
        }

        // notify the probe system
        probeSystem()->probeSecondary(iter);

        // verify and log loop convergence
        if (logLoopConvergence(log(), converged, iter, minIters, maxIters)) break;
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
                    // --- forced scattering ---
                    if (_config->forceScattering())
                    {
                        double Lthreshold = pp.luminosity() / _config->minWeightReduction();
                        int minScattEvents = _config->minScattEvents();
                        while (true)
                        {
                            // calculate segments and optical depths for the complete path
                            if (_config->explicitAbsorption())
                                mediumSystem()->setScatteringAndAbsorptionOpticalDepths(&pp);
                            else
                                mediumSystem()->setExtinctionOpticalDepths(&pp);

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
                    // --- non-forced scattering ---
                    else
                    {
                        while (true)
                        {
                            // advance the packet (without storing the radiation field)
                            // if the interaction point is outside of the path, terminate the packet
                            if (!simulateNonForcedPropagation(&pp)) break;

                            // if the packet's weight drops to zero, terminate it
                            if (pp.luminosity() <= 0) break;

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
                double lnExtEnd = -segment.tauExt();  // extinction factor and its logarithm at end of current segment
                double extEnd = exp(lnExtEnd);
                int m = segment.m();
                if (m >= 0)
                {
                    // use this flavor of the lnmean function to avoid recalculating the logarithm of the extinction
                    double extMean = SpecialFunctions::lnmean(extEnd, extBeg, lnExtEnd, lnExtBeg);
                    double Lds = luminosity * extMean * segment.ds();
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
            double lnExtEnd = -segment.tauExt();  // extinction factor and its logarithm at end of current segment
            double extEnd = exp(lnExtEnd);
            int m = segment.m();
            if (m >= 0)
            {
                double lambda = pp->perceivedWavelength(mediumSystem()->bulkVelocity(m),
                                                        _config->hubbleExpansionRate() * segment.s());
                int ell = _config->radiationFieldWLG()->bin(lambda);
                if (ell >= 0)
                {
                    // use this flavor of the lnmean function to avoid recalculating the logarithm of the extinction
                    double extMean = SpecialFunctions::lnmean(extEnd, extBeg, lnExtEnd, lnExtBeg);
                    double Lds = pp->perceivedLuminosity(lambda) * extMean * segment.ds();
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

    // adjust the photon packet weight with the escape fraction and, depending on the type of photon cycle,
    // with either the scattered fraction or the cumulative absorption optical depth at the interaction point
    if (_config->explicitAbsorption())
    {
        double tauAbs = pp->interactionOpticalDepth();
        pp->applyBias(-expm1(-taupath) * exp(-tauAbs));
    }
    else
    {
        double albedo = mediumSystem()->albedoForScattering(pp);
        pp->applyBias(-expm1(-taupath) * albedo);
    }

    // advance the photon packet position
    pp->propagate(pp->interactionDistance());
}

////////////////////////////////////////////////////////////////////

bool MonteCarloSimulation::simulateNonForcedPropagation(PhotonPacket* pp)
{
    // generate a random interaction optical depth
    double tauinteract = random()->expon();

    if (_config->explicitAbsorption())
    {
        // find the physical interaction point corresponding to this scattering optical depth
        // and calculate the absorption optical depth at the interaction point;
        // if the interaction point is outside of the path, terminate the photon packet
        if (!mediumSystem()->setInteractionPointUsingScatteringAndAbsorption(pp, tauinteract)) return false;

        // get the cumulative absorption optical depth at the interaction point
        double tauAbs = pp->interactionOpticalDepth();

        // adjust the photon packet weight by the corresponding extinction (or stimulation, for negative optical depth)
        pp->applyBias(exp(-tauAbs));
    }
    else
    {
        // find the physical interaction point corresponding to this optical depth
        // if the interaction point is outside of the path, terminate the photon packet
        if (!mediumSystem()->setInteractionPointUsingExtinction(pp, tauinteract)) return false;

        // calculate the albedo for the cell containing the interaction point
        double albedo = mediumSystem()->albedoForScattering(pp);

        // adjust the photon packet weight by the albedo
        pp->applyBias(albedo);
    }

    // advance the photon packet position
    pp->propagate(pp->interactionDistance());
    return true;
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

    // now do the actual peel-off
    if (_config->hasScatteringDispersion())
    {
        // if wavelengths may change, send a peel-off photon packet per medium component to each instrument
        int numMedia = wv.size();
        for (int h = 0; h != numMedia; ++h)
        {
            // skip media that don't scatter this photon packet
            if (wv[h] > 0.)
            {
                for (Instrument* instr : _instrumentSystem->instruments())
                {
                    if (!instr->isSameObserverAsPreceding())
                    {
                        // get the direction towards the instrument and (for polarization only) its Y-axis orientation
                        Direction bfkobs = instr->bfkobs(pp->position());
                        Direction bfky = _config->hasPolarization() ? instr->bfky(pp->position()) : Direction();

                        // calculate peel-off for the current component and launch the peel-off photon packet
                        mediumSystem()->peelOffScattering(h, wv[h], lambda, bfkobs, bfky, pp, ppp);
                    }

                    // have the peel-off photon packet detected
                    instr->detect(ppp);
                }
            }
        }
    }
    else
    {
        // if wavelengths cannot change, send a consolidated peel-off photon packet to each instrument
        for (Instrument* instr : _instrumentSystem->instruments())
        {
            if (!instr->isSameObserverAsPreceding())
            {
                // get the direction towards the instrument and (for polarization only) its Y-axis orientation
                Direction bfkobs = instr->bfkobs(pp->position());
                Direction bfky = _config->hasPolarization() ? instr->bfky(pp->position()) : Direction();

                // calculate peel-off for all medium components and launch the peel-off photon packet
                // (all media must either support polarization or not; combining these support levels is not allowed)
                mediumSystem()->peelOffScattering(wv, lambda, bfkobs, bfky, pp, ppp);
            }

            // have the peel-off photon packet detected
            instr->detect(ppp);
        }
    }
}

////////////////////////////////////////////////////////////////////
