/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustEmissionWavelengthGridProbe.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "InstrumentWavelengthGridProbe.hpp"
#include "MediumSystem.hpp"

////////////////////////////////////////////////////////////////////

void DustEmissionWavelengthGridProbe::probe()
{
    if (find<Configuration>()->hasDustEmission() && find<MediumSystem>()->hasDust())
    {
        InstrumentWavelengthGridProbe::writeWavelengthGrid(this, find<Configuration>()->dustEmissionWLG(),
                                                           itemName() + "_wavelengths",
                                                           "emission spectrum wavelengths");
    }
}

////////////////////////////////////////////////////////////////////
