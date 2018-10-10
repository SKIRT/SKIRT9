/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MonteCarloSimulation.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "MaterialMix.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PhotonPacket.hpp"
#include "ProcessManager.hpp"
#include "ShortArray.hpp"
#include "SpatialGrid.hpp"
#include "SpecialFunctions.hpp"
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

    // if there are no media, simply log the source model symmetry
    if (!_config->hasMedium())
    {
        log()->info("Model symmetry is " + std::to_string(sourceSystem()->dimension()) + "D");
    }

    // if there are media, compare the model symmetry to the grid symmetry
    else
    {
        int modelDimension = max(sourceSystem()->dimension(), mediumSystem()->dimension());
        int gridDimension = mediumSystem()->gridDimension();

        if (modelDimension == gridDimension)
        {
            log()->info("Model and grid symmetry: " + std::to_string(modelDimension) + "D");
        }
        else
        {
            log()->info("Model symmetry: " + std::to_string(modelDimension) + "D; "
                        "Spatial grid symmetry: " + std::to_string(gridDimension) + "D");
             if (modelDimension > gridDimension)
                 throw FATALERROR("This grid symmetry does not support the model symmetry");
             else
                 log()->warning("Selecting a grid with the model symmetry might be more efficient");
        }
    }
}

////////////////////////////////////////////////////////////////////

Configuration* MonteCarloSimulation::config() const
{
    return _config;
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

    log()->info("Launching " + StringUtils::toString(static_cast<double>(numTotal))
                + " " + _segment + " photon packages");
    log()->infoSetElapsed(numTotal);
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::logProgress(size_t numDone)
{
    // log message if the minimum time has elapsed
    log()->infoIfElapsed("Launched " + _segment + " photon packages: ", numDone);
}

////////////////////////////////////////////////////////////////////

namespace
{
    // maximum number of photon packets processed between two invocations of infoIfElapsed()
    const size_t logProgressChunkSize = 10000;
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

        // loop over the history indices, with interruptions for progress logging
        while (numIndices)
        {
            size_t currentChunkSize = min(logProgressChunkSize, numIndices);
            for (size_t historyIndex=firstIndex; historyIndex!=firstIndex+currentChunkSize; ++historyIndex)
            {
                // launch a photon packet
                sourceSystem()->launch(&pp, historyIndex);
                if (pp.luminosity()>0)
                {
                    peelOffEmission(&pp, &ppp);

                    // trace the packet through the media
                    double Lthreshold = pp.luminosity() / _config->minWeightReduction();
                    int minScattEvents = _config->minScattEvents();
                    if (_config->hasMedium()) while (true)
                    {
                        mediumSystem()->fillOpticalDepthInfo(&pp);
                        if (_config->hasRadiationField()) storeRadiationField(&pp);
                        simulatePropagation(&pp);
                        if (pp.luminosity()<=0 || (pp.luminosity()<=Lthreshold && pp.numScatt()>=minScattEvents)) break;
                        peelOffScattering(&pp,&ppp);
                        simulateScattering(&pp);
                    }
                }

            }

            // log progress
            logProgress(currentChunkSize);
            firstIndex += currentChunkSize;
            numIndices -= currentChunkSize;
        }
    }

#ifdef LOG_CHUNKS
    // time one of the threads for debugging purposes
    if (id==1) StopWatch::stop();
#endif
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::peelOffEmission(const PhotonPacket* pp, PhotonPacket* ppp)
{
    for (Instrument* instrument : _instrumentSystem->instruments())
    {
        ppp->launchEmissionPeelOff(pp, instrument->bfkobs(pp->position()));
        instrument->detect(ppp);
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::storeRadiationField(const PhotonPacket* pp)
{
    double extBeg = 1.;                         // extinction factor at begin of current segment
    for (const auto& segment : pp->segments())
    {
        double extEnd = exp(-segment.tau);      // extinction factor at end of current segment
        int m = segment.m;
        if (m >= 0)
        {
            double lambda = pp->perceivedWavelength(mediumSystem()->bulkVelocity(m));
            int ell = _config->radiationFieldWavelengthGrid()->bin(lambda);
            double Lds = pp->perceivedLuminosity(lambda) * SpecialFunctions::lnmean(extBeg, extEnd) * segment.ds;
            mediumSystem()->storeRadiationField(m, ell, Lds);
        }
        extBeg = extEnd;
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
    if (xi==0.)
    {
        tau = random()->exponCutoff(taupath);
    }
    else
    {
        tau = random()->uniform()<xi ? random()->uniform()*taupath : random()->exponCutoff(taupath);
        double p = -exp(-tau)/expm1(-taupath);
        double q = (1.0-xi)*p + xi/taupath;
        double weight = p/q;
        pp->applyBias(weight);
    }

    // determine the physical position of the interaction point
    pp->findInteractionPoint(tau);
    int m = pp->interactionCellIndex();
    if (m<0) throw FATALERROR("Cannot locate photon packet interaction point");

    // calculate the albedo for the cell containing the interaction point
    Vec bfv = mediumSystem()->bulkVelocity(m);
    double albedo = mediumSystem()->albedo(pp->perceivedWavelength(bfv), m);

    // adjust the weight by the scattered fraction
    pp->applyBias( -expm1(-taupath) * albedo );

    // advance the position
    pp->propagate(pp->interactionDistance());
}

////////////////////////////////////////////////////////////////////

namespace
{
    // This helper function returns the angle phi between the previous and current scattering planes
    // given the normal to the previous scattering plane and the current and new propagation directions
    // of the photon package. The function returns a zero angle if the light is unpolarized or when the
    // current scattering event is completely forward or backward.
    double angleBetweenScatteringPlanes(Direction np, Direction kc, Direction kn)
    {
        Vec nc = Vec::cross(kc,kn);
        nc /= nc.norm();
        double cosphi = Vec::dot(np,nc);
        double sinphi = Vec::dot(Vec::cross(np,nc), kc);
        double phi = atan2(sinphi,cosphi);
        if (std::isfinite(phi)) return phi;
        return 0.;
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::peelOffScattering(const PhotonPacket* pp, PhotonPacket* ppp)
{
    // get the cell hosting the scattering event
    int m = pp->interactionCellIndex();

    // get the bulk velocity of the material in that cell
    Vec bfv = mediumSystem()->bulkVelocity(m);

    // determine the perceived wavelength, Doppler-shifted from the incoming wavelength according to the bulk velocity
    double lambda = pp->perceivedWavelength(bfv);

    // determine the weighting factor for each medium component as its scattering opacity (n * sigma_sca)
    int numMedia = mediumSystem()->numMedia();
    ShortArray<8> wv(numMedia);
    if (numMedia==1)
    {
        wv[0] = 1.;
    }
    else
    {
        double sum = 0.;
        for (int h=0; h!=numMedia; ++h)
        {
            wv[h] = mediumSystem()->opacitySca(lambda, m, h);
            sum += wv[h];
        }
        if (sum<=0) return; // abort peel-off if none of the media scatters
        for (int h=0; h!=numMedia; ++h) wv[h] /= sum;
    }

    // now do the actual peel-off for each instrument
    for (Instrument* instr : _instrumentSystem->instruments())
    {
        // get the instrument direction
        Direction bfkobs = instr->bfkobs(pp->position());

        // calculate the weighted sum of the effects on the Stokes vector for all media
        double I = 0., Q = 0., U = 0., V = 0.;
        for (int h=0; h!=numMedia; ++h)
        {
            // use the appropriate algorithm for each mix
            // (all mixes must either support polarization or not; combining these support levels is not allowed)
            auto mix = mediumSystem()->mix(m,h);
            switch (mix->scatteringMode())
            {
            case MaterialMix::ScatteringMode::HenyeyGreenstein:
                {
                    // calculate the value of the Henyey-Greenstein phase function
                    double costheta = Vec::dot(pp->direction(), bfkobs);
                    double g = mix->asymmpar(lambda);
                    double t = 1.0+g*g-2*g*costheta;
                    double value = (1.0-g)*(1.0+g)/sqrt(t*t*t);

                    // accumulate the weighted sum in the intensity (there is no support for polarization in this case)
                    I += wv[h] * value;
                    break;
                }
            case MaterialMix::ScatteringMode::MaterialPhaseFunction:
                {
                    // calculate the value of the material-specific phase function
                    double costheta = Vec::dot(pp->direction(), bfkobs);
                    double value = mix->phaseFunctionValueForCosine(lambda, costheta);

                    // accumulate the weighted sum in the intensity (there is no support for polarization in this case)
                    I += wv[h] * value;
                    break;
                }
            case MaterialMix::ScatteringMode::SphericalPolarization:
                {
                    // calculate the value of the material-specific phase function
                    double theta = acos(Vec::dot(pp->direction(),bfkobs));
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
        instr->detect(ppp);
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::simulateScattering(PhotonPacket* pp)
{
    // locate the cell hosting the scattering event
    int m = pp->interactionCellIndex();

    // get the bulk velocity of the material in that cell
    Vec bfv = mediumSystem()->bulkVelocity(m);

    // determine the perceived wavelength, Doppler-shifted from the incoming wavelength according to the bulk velocity
    double lambda = pp->perceivedWavelength(bfv);

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
                double f = ((1.0-g)*(1.0+g))/(1.0-g+2.0*g*random()->uniform());
                double costheta = (1.0+g*g-f*f)/(2.0*g);
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
            Vec newdir = pp->direction()*cos(theta) + Vec::cross(pp->normal(), pp->direction())*sin(theta);

            // normalize the new direction to prevent degradation
            bfknew = Direction(newdir/newdir.norm());
            break;
        }
    }
    pp->scatter(bfknew, bfv);
}

////////////////////////////////////////////////////////////////////
