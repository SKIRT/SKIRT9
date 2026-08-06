/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpatialGridWhenFormProbe.hpp"

////////////////////////////////////////////////////////////////////

Probe::When SpatialGridWhenFormProbe::when() const
{
    switch (probeAfter())
    {
        case ProbeAfter::Setup: return When::Setup;
        case ProbeAfter::Run: return When::Run;
        case ProbeAfter::Primary: return When::Primary;
        case ProbeAfter::Secondary: return When::Secondary;
    }
    return When::Setup;
}

////////////////////////////////////////////////////////////////////
