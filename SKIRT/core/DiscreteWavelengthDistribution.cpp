/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DiscreteWavelengthDistribution.hpp"
#include "FatalError.hpp"
#include "Random.hpp"
#include "SourceWavelengthRangeInterface.hpp"

//////////////////////////////////////////////////////////////////////

void DiscreteWavelengthDistribution::setupSelfAfter()
{
    WavelengthDistribution::setupSelfAfter();

    // get the source wavelength range
    Range range = interface<SourceWavelengthRangeInterface>()->wavelengthRange();

    // determine the wavelengths that fall inside the source range
    int numBins = wavelengthGrid()->numBins();
    for (_beginWavelengthIndex = 0; _beginWavelengthIndex != numBins; ++_beginWavelengthIndex)
    {
        if (range.contains(wavelengthGrid()->wavelength(_beginWavelengthIndex))) break;
    }
    for (_endWavelengthIndex = _beginWavelengthIndex; _endWavelengthIndex != numBins; ++_endWavelengthIndex)
    {
        if (!range.contains(wavelengthGrid()->wavelength(_endWavelengthIndex))) break;
    }
    _numWavelengths = _endWavelengthIndex - _beginWavelengthIndex;

    // verify that there is at least one wavelength and that the smallest one is positive
    if (_numWavelengths <= 0) throw FATALERROR("None of the grid wavelengths are in the source wavelength range");
}

//////////////////////////////////////////////////////////////////////

double DiscreteWavelengthDistribution::probability(double wavelength) const
{
    // get the wavelength index and verify that it is within range
    int ell = wavelengthGrid()->bin(wavelength);
    if (ell < _beginWavelengthIndex || ell >= _endWavelengthIndex) return 0.;

    // calculate the probability for wavelengths within this bin
    return 1. / (_numWavelengths * wavelengthGrid()->effectiveWidth(ell));
}

//////////////////////////////////////////////////////////////////////

double DiscreteWavelengthDistribution::generateWavelength() const
{
    int ell = _beginWavelengthIndex + static_cast<int>(random()->uniform() * _numWavelengths);
    return wavelengthGrid()->wavelength(ell);
}

//////////////////////////////////////////////////////////////////////
