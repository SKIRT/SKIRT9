/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ListWavelengthGrid.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

void ListWavelengthGrid::setupSelfBefore()
{
    DisjointWavelengthGrid::setupSelfBefore();

    // set the wavelength grid from the list of property values
    if (_relativeHalfWidth)
        setWavelengthBins(NR::array(_wavelengths), _relativeHalfWidth);
    else
        setWavelengthRange(NR::array(_wavelengths), _log);
}

//////////////////////////////////////////////////////////////////////
