/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LyaSEDDecorator.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void LyaSEDDecorator::setupSelfAfter()
{
    SED::setupSelfAfter();

    // determine the ionized luminosity fraction
    Range totalRange = normalizationWavelengthRange();
    if (totalRange.min() >= Constants::lambdaIon())
        throw FATALERROR("Source wavelength range must include ionizing radiation");
    _ionizingFraction = _sedOriginal->integratedLuminosity(Range(totalRange.min(), Constants::lambdaIon()))
                        / _sedOriginal->integratedLuminosity(totalRange);
}

//////////////////////////////////////////////////////////////////////

Range LyaSEDDecorator::intrinsicWavelengthRange() const
{
    Range range = sedOriginal()->intrinsicWavelengthRange();
    range.extend(sedLymanAlpha()->intrinsicWavelengthRange());
    return range;
}

//////////////////////////////////////////////////////////////////////

double LyaSEDDecorator::specificLuminosity(double wavelength) const
{
    if (wavelength <= Constants::lambdaIon())
        return (1. - conversionFraction()) * sedOriginal()->specificLuminosity(wavelength);
    else
        return sedOriginal()->specificLuminosity(wavelength)
               + _ionizingFraction * conversionFraction() * sedLymanAlpha()->specificLuminosity(wavelength);
}

//////////////////////////////////////////////////////////////////////

void LyaSEDDecorator::specificLuminosityArray(Array& lambdav, Array& pv, const Range& wavelengthRange) const
{
    throw FATALERROR("not yet implemented");
}

//////////////////////////////////////////////////////////////////////

double LyaSEDDecorator::integratedLuminosity(const Range& wavelengthRange) const
{
    double luminosity = 0;

    Range ionizingRange(std::numeric_limits<double>::denorm_min(), Constants::lambdaIon());
    ionizingRange.intersect(wavelengthRange);
    if (!ionizingRange.empty())
    {
        luminosity += (1. - conversionFraction()) * sedOriginal()->integratedLuminosity(ionizingRange);
    }

    Range nonIonizingRange(Constants::lambdaIon(), std::numeric_limits<double>::max());
    nonIonizingRange.intersect(wavelengthRange);
    if (!nonIonizingRange.empty())
    {
        luminosity +=
            sedOriginal()->integratedLuminosity(nonIonizingRange)
            + _ionizingFraction * conversionFraction() * sedLymanAlpha()->integratedLuminosity(nonIonizingRange);
    }

    return luminosity;
}

//////////////////////////////////////////////////////////////////////

double LyaSEDDecorator::generateWavelength() const
{
    double wavelength = sedOriginal()->generateWavelength();
    if (wavelength <= Constants::lambdaIon())
    {
        if (random()->uniform() <= conversionFraction())
        {
            wavelength = sedLymanAlpha()->generateWavelength();
        }
    }
    return wavelength;
}

//////////////////////////////////////////////////////////////////////
