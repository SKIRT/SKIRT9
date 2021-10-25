/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "EmittingGasMix.hpp"
#include "DisjointWavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

Range EmittingGasMix::wavelengthRange() const
{
    Range lineRange, contRange;
    if (hasLineEmission())
    {
        Array centers = lineEmissionCenters();
        lineRange.set(centers.min()*0.999, centers.max()*1.001);
    }
    if (hasContinuumEmission())
    {
        contRange = continuumEmissionWLG()->wavelengthRange();
    }

    if (lineRange.empty()) return contRange;
    if (contRange.empty()) return lineRange;
    return lineRange.extend(contRange);
}

////////////////////////////////////////////////////////////////////
