/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MonteCarloSimulation.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "MaterialMix.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PhotonPacket.hpp"
#include "ProcessManager.hpp"
#include "SecondarySourceSystem.hpp"
#include "ShortArray.hpp"
#include "SpatialGrid.hpp"
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

        // opacity iteration. TODO: enforce requirements in Configuration once we make up our minds
        bool withSecondary = false;
        if (mediumSystem()->hasGas() && _config->hasRadiationField())
        {
            // Iteration without secondary emission (big changes)
            runSelfConsistentOpacityPhase(false);
            // Iteration with secondary emission (smaller changes)
            if (withSecondary) runSelfConsistentOpacityPhase(true);
        }

        // primary emission segment
        runPrimaryEmission();

        // dust self-absorption iteration segments (still needs to run if opacity iteration did not
        // treat secondary emission)
        if (_config->hasDustSelfAbsorption() && !withSecondary) runDustSelfAbsorptionPhase();

        // secondary emission segment
        if (_config->hasSecondaryEmission()) runSecondaryEmission();

        _mediumSystem->gasTest();
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
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runSelfConsistentOpacityPhase(bool)
{
    TimeLogger logger(log(), "the self consistent opacity phase");

    size_t Npp = _config->numPrimaryPackets();
    if (!Npp)
    {
        log()->warning("Skipping primary emission because no photon packets were requested");
        return;
    }
    else if (!sourceSystem()->luminosity())
    {
        log()->warning("Skipping primary emission because the total luminosity of primary sources is zero");
        return;
    }

    auto parallel = find<ParallelFactory>()->parallelDistributed();

    int maxIters = 10;

    // initialize some convergence criteria here

    for (int iter = 1; iter <= maxIters; iter++)
    {
        string segment = "self consistent opacity iteration " + std::to_string(iter);
        {
            TimeLogger logger(log(), segment);

            mediumSystem()->clearRadiationField(true);

            initProgress(segment, Npp);
            sourceSystem()->prepareForLaunch(Npp);
            parallel->call(Npp, [this](size_t i, size_t n) { performLifeCycle(i, n, true, false, true); });
            instrumentSystem()->flush();
            wait(segment);
            mediumSystem()->communicateRadiationField(true);

            // initialize or update the gas opacities using the stored radiation field
            mediumSystem()->updateGas();
        }
        // check convergence here
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runDustSelfAbsorptionPhase()
{
    TimeLogger logger(log(), "the dust self-absorption phase");

    // get number of photons; return if zero
    size_t Npp = _config->numIterationPackets();
    if (!Npp)
    {
        log()->warning("Skipping dust self-absorption phase because no photon packets were requested");
        return;
    }

    // get the parallel engine
    auto parallel = find<ParallelFactory>()->parallelDistributed();

    // get the parameters controlling the self-absorption iteration
    int minIters = _config->minIterations();
    int maxIters = _config->maxIterations();
    double fractionOfPrimary = _config->maxFractionOfPrimary();
    double fractionOfPrevious = _config->maxFractionOfPrevious();

    // initialize the total absorbed luminosity in the previous iteration
    double prevLabsdust = 0.;

    // iterate over the maximum number of iterations; the loop body returns from the function
    // when convergence is reached after the minimum number of iterations have been completed
    for (int iter = 1; iter <= maxIters; iter++)
    {
        string segment = "dust self-absorption iteration " + std::to_string(iter);
        {
            TimeLogger logger(log(), segment);

            // clear the secondary radiation field
            mediumSystem()->clearRadiationField(false);

            // prepare the source system; terminate if the dust has zero luminosity (which should never happen)
            if (!_secondarySourceSystem->prepareForLaunch(Npp))
            {
                log()->warning("Terminating dust self-absorption phase because the total dust luminosity is zero");
                return;
            }

            // launch photon packets
            initProgress(segment, Npp);
            parallel->call(Npp, [this](size_t i, size_t n) { performLifeCycle(i, n, false, false, true); });
            instrumentSystem()->flush();

            // wait for all processes to finish and synchronize the radiation field
            wait(segment);
            mediumSystem()->communicateRadiationField(false);
        }

        // determine and log the total absorbed luminosity
        double Labsprim = mediumSystem()->totalAbsorbedLuminosity(true, MaterialMix::MaterialType::Dust);
        double Labsdust = mediumSystem()->totalAbsorbedLuminosity(false, MaterialMix::MaterialType::Dust);
        log()->info("The total dust-absorbed primary luminosity is "
                    + StringUtils::toString(units()->obolluminosity(Labsprim), 'g') + " " + units()->ubolluminosity());
        log()->info("The total dust-absorbed dust luminosity in iteration " + std::to_string(iter) + " is "
                    + StringUtils::toString(units()->obolluminosity(Labsdust), 'g') + " " + units()->ubolluminosity());

        // log the current performance and corresponding convergence criteria
        if (Labsprim > 0. && Labsdust > 0.)
        {
            if (iter == 1)
            {
                log()->info("--> absorbed dust luminosity is "
                            + StringUtils::toString(Labsdust / Labsprim * 100., 'f', 2)
                            + "% of absorbed stellar luminosity (convergence criterion is "
                            + StringUtils::toString(fractionOfPrimary * 100., 'f', 2) + "%)");
            }
            else
            {
                log()->info("--> absorbed dust luminosity changed by "
                            + StringUtils::toString(abs((Labsdust - prevLabsdust) / Labsdust) * 100., 'f', 2)
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
            // - the absorbed stellar luminosity is zero
            // - the absorbed dust luminosity is zero
            // - the absorbed dust luminosity is less than a given fraction of the absorbed stellar luminosity
            // - the absorbed dust luminosity has changed by less than a given fraction compared to the previous iter
            if (Labsprim <= 0. || Labsdust <= 0. || Labsdust / Labsprim < fractionOfPrimary
                || abs((Labsdust - prevLabsdust) / Labsdust) < fractionOfPrevious)
            {
                log()->info("Convergence reached after " + std::to_string(iter) + " iterations");
                return;  // end the iteration by returning from the function
            }
            else
            {
                log()->info("Convergence not yet reached after " + std::to_string(iter) + " iterations");
            }
        }
        prevLabsdust = Labsdust;
    }

    // if the loop runs out, convergence was not reached even after the maximum number of iterations
    log()->error("Convergence not yet reached after " + std::to_string(maxIters) + " iterations");
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runSecondaryEmission()
{
    string segment = "secondary emission";
    TimeLogger logger(log(), segment);

    // determine whether we need to store the radiation field during secondary emission
    bool storeRF = _config->storeEmissionRadiationField();

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
    const size_t logProgressChunkSize = 10000;
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
                    double Lthreshold = pp.luminosity() / _config->minWeightReduction();
                    int minScattEvents = _config->minScattEvents();
                    while (true)
                    {
                        mediumSystem()->opticalDepth(&pp);
                        if (store) storeRadiationField(&pp);
                        simulatePropagation(&pp);
                        if (pp.luminosity() <= 0 || (pp.luminosity() <= Lthreshold && pp.numScatt() >= minScattEvents))
                            break;
                        if (peel) peelOffScattering(&pp, &ppp);
                        simulateScattering(&pp);
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
    if (!_config->hasMovingMedia())
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
                double lambda = pp->perceivedWavelength(mediumSystem()->bulkVelocity(m));
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

void MonteCarloSimulation::simulatePropagation(PhotonPacket* pp)
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
    int m = pp->interactionCellIndex();
    if (m < 0) throw FATALERROR("Cannot locate photon packet interaction point");

    // calculate the albedo for the cell containing the interaction point
    // use a faster version in case there are no kinematics
    double albedo;
    if (!_config->hasMovingMedia())
    {
        albedo = mediumSystem()->albedo(pp->wavelength(), m);
    }
    else
    {
        Vec bfv = mediumSystem()->bulkVelocity(m);
        albedo = mediumSystem()->albedo(pp->perceivedWavelength(bfv), m);
    }

    // adjust the weight by the scattered fraction
    pp->applyBias(-expm1(-taupath) * albedo);

    // advance the position
    pp->propagate(pp->interactionDistance());
}

////////////////////////////////////////////////////////////////////

namespace
{
    // This helper function returns the angle phi between the previous and current scattering planes
    // given the normal to the previous scattering plane and the current and new propagation directions
    // of the photon packet. The function returns a zero angle if the light is unpolarized or when the
    // current scattering event is completely forward or backward.
    double angleBetweenScatteringPlanes(Direction np, Direction kc, Direction kn)
    {
        Vec nc = Vec::cross(kc, kn);
        nc /= nc.norm();
        double cosphi = Vec::dot(np, nc);
        double sinphi = Vec::dot(Vec::cross(np, nc), kc);
        double phi = atan2(sinphi, cosphi);
        if (std::isfinite(phi)) return phi;
        return 0.;
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::peelOffScattering(const PhotonPacket* pp, PhotonPacket* ppp)
{
    // get the cell hosting the scattering event
    int m = pp->interactionCellIndex();

    // get the bulk velocity of the material in that cell and determine the perceived wavelength,
    // Doppler-shifted from the incoming wavelength according to the bulk velocity;
    // use a faster version in case there are no kinematics
    Vec bfv;
    double lambda;
    if (!_config->hasMovingMedia())
    {
        lambda = pp->wavelength();
    }
    else
    {
        bfv = mediumSystem()->bulkVelocity(m);
        lambda = pp->perceivedWavelength(bfv);
    }

    // determine the weighting factor for each medium component as its scattering opacity (n * sigma_sca)
    int numMedia = mediumSystem()->numMedia();
    ShortArray<8> wv(numMedia);
    if (numMedia == 1)
    {
        wv[0] = 1.;
    }
    else
    {
        double sum = 0.;
        for (int h = 0; h != numMedia; ++h)
        {
            wv[h] = mediumSystem()->opacitySca(lambda, m, h);
            sum += wv[h];
        }
        if (sum <= 0) return;  // abort peel-off if none of the media scatters
        for (int h = 0; h != numMedia; ++h) wv[h] /= sum;
    }

    // now do the actual peel-off for each instrument
    for (Instrument* instr : _instrumentSystem->instruments())
    {
        if (!instr->isSameObserverAsPreceding())
        {
            // get the instrument direction
            Direction bfkobs = instr->bfkobs(pp->position());

            // calculate the weighted sum of the effects on the Stokes vector for all media
            double I = 0., Q = 0., U = 0., V = 0.;
            for (int h = 0; h != numMedia; ++h)
            {
                // use the appropriate algorithm for each mix
                // (all mixes must either support polarization or not; combining these support levels is not allowed)
                auto mix = mediumSystem()->mix(m, h);
                switch (mix->scatteringMode())
                {
                    case MaterialMix::ScatteringMode::HenyeyGreenstein:
                    {
                        // calculate the value of the Henyey-Greenstein phase function
                        double costheta = Vec::dot(pp->direction(), bfkobs);
                        double g = mix->asymmpar(lambda);
                        double t = 1.0 + g * g - 2 * g * costheta;
                        double value = (1.0 - g) * (1.0 + g) / sqrt(t * t * t);

                        // accumulate the weighted sum in the intensity (no support for polarization in this case)
                        I += wv[h] * value;
                        break;
                    }
                    case MaterialMix::ScatteringMode::MaterialPhaseFunction:
                    {
                        // calculate the value of the material-specific phase function
                        double costheta = Vec::dot(pp->direction(), bfkobs);
                        double value = mix->phaseFunctionValueForCosine(lambda, costheta);

                        // accumulate the weighted sum in the intensity (no support for polarization in this case)
                        I += wv[h] * value;
                        break;
                    }
                    case MaterialMix::ScatteringMode::SphericalPolarization:
                    case MaterialMix::ScatteringMode::SpheroidalPolarization:
                    {
                        // calculate the value of the material-specific phase function
                        double theta = acos(Vec::dot(pp->direction(), bfkobs));
                        double phi = angleBetweenScatteringPlanes(pp->normal(), pp->direction(), bfkobs);
                        double value = mix->phaseFunctionValue(lambda, theta, phi, pp);

                        // copy the polarization state so we can change it without affecting the incoming photon packet
                        StokesVector sv = *pp;

                        // rotate the Stokes vector reference direction into the scattering plane
                        sv.rotateIntoPlane(pp->direction(), bfkobs);

                        // apply the Mueller matrix
                        mix->applyMueller(lambda, theta, &sv);

                        // rotate the Stokes vector reference direction parallel to the instrument frame y-axis
                        // it is given bfkobs because the photon is at this point aimed towards the observer
                        sv.rotateIntoPlane(bfkobs, instr->bfky(pp->position()));

                        // acumulate the weighted sum of all Stokes components to support polarization
                        double w = wv[h] * value;
                        I += w * sv.stokesI();
                        Q += w * sv.stokesQ();
                        U += w * sv.stokesU();
                        V += w * sv.stokesV();
                        break;
                    }
                }
            }

            // pass the result to the peel-off photon packet and have it detected
            ppp->launchScatteringPeelOff(pp, bfkobs, bfv, I);
            if (_config->hasPolarization()) ppp->setPolarized(I, Q, U, V, pp->normal());
        }
        instr->detect(ppp);
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::simulateScattering(PhotonPacket* pp)
{
    // locate the cell hosting the scattering event
    int m = pp->interactionCellIndex();

    // get the bulk velocity of the material in that cell and determine the perceived wavelength,
    // Doppler-shifted from the incoming wavelength according to the bulk velocity;
    // use a faster version in case there are no kinematics
    Vec bfv;
    double lambda;
    if (!_config->hasMovingMedia())
    {
        lambda = pp->wavelength();
    }
    else
    {
        bfv = mediumSystem()->bulkVelocity(m);
        lambda = pp->perceivedWavelength(bfv);
    }

    // randomly select a material mix; the probability of each component is weighted by the scattering opacity
    auto mix = mediumSystem()->randomMixForScattering(random(), lambda, m);

    // now perform the scattering using this material mix
    //   - determine the new propagation direction
    //   - if supported, update the polarization state of the photon packet along the way
    Direction bfknew;
    switch (mix->scatteringMode())
    {
        case MaterialMix::ScatteringMode::HenyeyGreenstein:
        {
            // sample a scattering angle from the Henyey-Greenstein phase function
            // handle isotropic scattering separately because the HG sampling procedure breaks down in this case
            double g = mix->asymmpar(lambda);
            if (fabs(g) < 1e-6)
            {
                bfknew = random()->direction();
            }
            else
            {
                double f = ((1.0 - g) * (1.0 + g)) / (1.0 - g + 2.0 * g * random()->uniform());
                double costheta = (1.0 + g * g - f * f) / (2.0 * g);
                bfknew = random()->direction(pp->direction(), costheta);
            }
            break;
        }
        case MaterialMix::ScatteringMode::MaterialPhaseFunction:
        {
            // sample a scattering angle from the material-specific phase function
            double costheta = mix->generateCosineFromPhaseFunction(lambda);
            bfknew = random()->direction(pp->direction(), costheta);
            break;
        }
        case MaterialMix::ScatteringMode::SphericalPolarization:
        case MaterialMix::ScatteringMode::SpheroidalPolarization:
        {
            // sample the angles between the previous and new direction from the material-specific phase function,
            // given the incoming polarization state
            double theta, phi;
            std::tie(theta, phi) = mix->generateAnglesFromPhaseFunction(lambda, pp);

            // rotate the Stokes vector (and the scattering plane) of the photon packet
            pp->rotateStokes(phi, pp->direction());

            // apply Mueller matrix to the Stokes vector of the photon packet
            mix->applyMueller(lambda, theta, pp);

            // rotate the propagation direction in the scattering plane
            Vec newdir = pp->direction() * cos(theta) + Vec::cross(pp->normal(), pp->direction()) * sin(theta);

            // normalize the new direction to prevent degradation
            bfknew = Direction(newdir / newdir.norm());
            break;
        }
    }
    pp->scatter(bfknew, bfv);
}

////////////////////////////////////////////////////////////////////
