/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SecondarySourceSystem.hpp"
#include "Configuration.hpp"
#include "EmittingGasMix.hpp"
#include "FatalError.hpp"
#include "MediumSystem.hpp"
#include "PhotonPacket.hpp"
#include "ProbePhotonPacketInterface.hpp"
#include "SecondaryDustSource.hpp"

////////////////////////////////////////////////////////////////////

SecondarySourceSystem::SecondarySourceSystem(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

void SecondarySourceSystem::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    auto config = find<Configuration>();
    auto ms = find<MediumSystem>();

    // add aggregate dust source if applicable
    if (config->hasDustEmission() && ms->hasDust())
    {
        _sources.push_back(new SecondaryDustSource(this));
    }

    // add gas sources if applicable
    if (config->hasGasEmission())
    {
        for (int h : ms->gasMediumIndices())
        {
            auto emittingMix = dynamic_cast<const EmittingGasMix*>(ms->mix(0, h));
            if (emittingMix)
            {
                //if (emittingMix->hasContinuumEmission()) _sources.push_back(new SecondaryDustSource(this));
                //if (emittingMix->hasLineEmission()) _sources.push_back(new SecondaryDustSource(this));
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

void SecondarySourceSystem::installLaunchCallBack(ProbePhotonPacketInterface* callback)
{
    if (_callback) throw FATALERROR("Cannot install more than one photon packet launch probe");
    _callback = callback;
}

////////////////////////////////////////////////////////////////////

bool SecondarySourceSystem::prepareForLaunch(size_t numPackets)
{
    // calculate the total luminosity over all sources; if it is zero, report failure
    double L = 0.;
    for (auto source : _sources) L += source->prepareLuminosities();
    if (L <= 0.) return false;

    // tell the sources to prepare their packet mappings
    for (auto source : _sources) source->preparePacketMap(0, numPackets);

    // report success
    return true;
}

////////////////////////////////////////////////////////////////////

void SecondarySourceSystem::launch(PhotonPacket* pp, size_t historyIndex) const
{
    // select the source from which to launch based on the history index of this photon packet
    int s = 0;

    // tell the source to launch a packet
    _sources[s]->launch(pp, historyIndex);

    // add origin info (the index reflects the secondary source, not the medium component)
    pp->setSecondaryOrigin(s);

    // invoke launch call-back if installed
    if (_callback) _callback->probePhotonPacket(pp);
}

////////////////////////////////////////////////////////////////////
