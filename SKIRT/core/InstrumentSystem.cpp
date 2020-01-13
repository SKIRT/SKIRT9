/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "InstrumentSystem.hpp"

////////////////////////////////////////////////////////////////////

void InstrumentSystem::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    Instrument* preceding = nullptr;
    for (Instrument* instrument : _instruments)
    {
        if (!preceding)
            preceding = instrument;
        else
            instrument->determineSameObserverAsPreceding(preceding);
    }
}

////////////////////////////////////////////////////////////////////

void InstrumentSystem::flush()
{
    for (Instrument* instrument : _instruments) instrument->flush();
}

////////////////////////////////////////////////////////////////////

void InstrumentSystem::write()
{
    for (Instrument* instrument : _instruments) instrument->write();
}

////////////////////////////////////////////////////////////////////
