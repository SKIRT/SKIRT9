/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "OligoWavelengthDistribution.hpp"
#include "NR.hpp"
#include "OligoWavelengthGrid.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

OligoWavelengthDistribution::OligoWavelengthDistribution(OligoWavelengthGrid* wavelengthGrid)
    // initialize data members
    : _wavelengths(wavelengthGrid->lambdav()),
      _probability(1. / wavelengthGrid->numBins() / wavelengthGrid->effectiveWidth(0))
{
    // hook ourselves into the simulation hierarchy
    wavelengthGrid->addChild(this);

    // perform further setup (perhaps in one of the base classes)
    setup();
}

////////////////////////////////////////////////////////////////////

double OligoWavelengthDistribution::probability(double /*wavelength*/) const
{
    return _probability;
}

////////////////////////////////////////////////////////////////////

double OligoWavelengthDistribution::generateWavelength() const
{
    size_t index = static_cast<size_t>(random()->uniform() * _wavelengths.size());
    return _wavelengths[index];
}

////////////////////////////////////////////////////////////////////
