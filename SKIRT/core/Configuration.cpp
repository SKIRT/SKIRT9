/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Configuration.hpp"
#include "ExtinctionOnlyMode.hpp"
#include "FatalError.hpp"
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


    // retrieve options from NoMediaMode
    {
        auto mode = find<NoMediaMode>(false);
        if (!mode) throw FATALERROR("Cannot locate a SimulationMode object in the simulation hierarchy");
        _numPrimaryPackets = mode->numPackets() * mode->primaryPacketsMultiplier();
    }

    // retrieve extra options from ExtinctionOnlyMode
    {
        auto mode = find<ExtinctionOnlyMode>(false);
        if (mode)
        {
            _minWeightReduction = mode->minWeightReduction();
            _minScattEvents = mode->minScattEvents();
            _pathLengthBias = mode->pathLengthBias();
        }
    }

    // TODO: verify configuration for the above two modes
    // TODO: check for polarization
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
