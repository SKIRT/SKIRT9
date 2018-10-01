/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "OligoWavelengthGrid.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

OligoWavelengthGrid::OligoWavelengthGrid(SimulationItem* parent, const vector<double>& wavelengths)
{
    parent->addChild(this);
    _wavelengths = wavelengths;
    setup();
}

////////////////////////////////////////////////////////////////////

void OligoWavelengthGrid::setupSelfBefore()
{
    DisjointWavelengthGrid::setupSelfBefore();

    // set the wavelength grid from the list of property values; use the same absolute width for all bins
    setWavelengthBins(NR::array(_wavelengths), 1e-3, true);
}

////////////////////////////////////////////////////////////////////
