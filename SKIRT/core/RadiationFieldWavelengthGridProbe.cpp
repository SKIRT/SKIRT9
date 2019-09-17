/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "RadiationFieldWavelengthGridProbe.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "MediumSystem.hpp"
#include "SpatialGrid.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"
#include "InstrumentWavelengthGridProbe.hpp"

////////////////////////////////////////////////////////////////////

void RadiationFieldWavelengthGridProbe::probeSetup()
{
    if (find<Configuration>()->hasRadiationField())
    {
        InstrumentWavelengthGridProbe::writeWavelengthGrid(this, find<Configuration>()->radiationFieldWLG(),
                                                 itemName() + "_wavelengths", "wavelengths for mean intensity");
    }
}

////////////////////////////////////////////////////////////////////
