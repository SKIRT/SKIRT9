/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Configuration.hpp"
#include "ExtinctionOnlyMode.hpp"
#include "FatalError.hpp"
#include "MaterialMix.hpp"
#include "MediumSystem.hpp"
#include "NoMediaMode.hpp"

////////////////////////////////////////////////////////////////////

Configuration::Configuration(SimulationItem* parent)
{
    parent->addChild(this);
}

////////////////////////////////////////////////////////////////////

void Configuration::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // Impementation note: this function is NOT allowed to perform setup on any simulation item in the hierarchy;
    //                     in other words, always use find<XXX>(false) and check for a nullptr result.

    // keep track of configuration requirements that need to be verified
    bool mustHaveMedia = false;

    // retrieve options from NoMediaMode
    {
        auto mode = find<NoMediaMode>(false);
        if (!mode) throw FATALERROR("Cannot locate a SimulationMode object in the simulation hierarchy");
        _numPrimaryPackets = mode->numPackets() * mode->primaryPacketsMultiplier();
    }

    // retrieve extra options from ExtinctionOnlyMode, if present
    {
        auto mode = find<ExtinctionOnlyMode>(false);
        if (mode)
        {
            mustHaveMedia = true;
            _minWeightReduction = mode->minWeightReduction();
            _minScattEvents = mode->minScattEvents();
            _pathLengthBias = mode->pathLengthBias();
            _numDensitySamples = mode->numDensitySamples();
        }
    }

    // determine the number of media in the simulation hierarchy
    int numMedia = 0;
    auto ms = find<MediumSystem>(false);
    if (ms) numMedia = ms->numMedia();  // may be zero

    // verify this with the requirements set by the simulation mode
    if (!mustHaveMedia && numMedia)
        throw FATALERROR("This simulation mode does not allow media to be configured");
    if (mustHaveMedia && !numMedia)
        throw FATALERROR("This simulation mode requires at least one medium to be configured");

    // check for polarization
    if (numMedia)
    {
        int numPolarization = 0;
        for (auto medium : ms->media()) if (medium->mix()->hasPolarization()) numPolarization++;
        if (numPolarization!=0 && numPolarization!=numMedia)
            throw FATALERROR("All media must consistenly support polarization, or not support polarization");
        _hasPolarization = numPolarization!=0;
    }

    // TODO: retrieve options from and verify configuration for other modes

    // in case emulation mode has been set before our setup() was called, perform the emulation overrides again
    if (emulationMode()) setEmulationMode();
}

////////////////////////////////////////////////////////////////////

void Configuration::setEmulationMode()
{
    _emulationMode = true;
    _numPrimaryPackets = 0.;
    // TODO: set other number of packet variables to zero
    // TODO: force the number of state iterations to one, if applicable
}

////////////////////////////////////////////////////////////////////
