/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ListBorderWavelengthGrid.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

void ListBorderWavelengthGrid::setupSelfBefore()
{
    DisjointWavelengthGrid::setupSelfBefore();

    // set the wavelength grid
    Array waves = NR::array(_wavelengths);
    switch (_characteristic)
    {
        case Characteristic::Linear: setWavelengthBorders(waves, false); break;
        case Characteristic::Logarithmic: setWavelengthBorders(waves, true); break;
        case Characteristic::Specified: setWavelengthSegments(waves); break;
    }
}

//////////////////////////////////////////////////////////////////////
