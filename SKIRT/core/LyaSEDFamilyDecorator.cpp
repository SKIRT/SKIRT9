/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LyaSEDFamilyDecorator.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

//////////////////////////////////////////////////////////////////////

void LyaSEDFamilyDecorator::setupSelfAfter()
{
    SEDFamily::setupSelfAfter();

    // determine the total and ionizing wavelength ranges
    Range range = find<Configuration>()->sourceWavelengthRange();
    range.intersect(sedFamilyOriginal()->intrinsicWavelengthRange());
    if (range.empty()) throw FATALERROR("Source wavelength range must overlap intrinsic SED family wavelength range");
    if (range.min() >= Constants::lambdaIon())
        throw FATALERROR("Source wavelength range must include ionizing radiation");
    _ionizingRange.set(range.min(), Constants::lambdaIon());
}

//////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> LyaSEDFamilyDecorator::parameterInfo() const
{
    return sedFamilyOriginal()->parameterInfo();
}

//////////////////////////////////////////////////////////////////////

Range LyaSEDFamilyDecorator::intrinsicWavelengthRange() const
{
    Range range = sedFamilyOriginal()->intrinsicWavelengthRange();
    range.extend(sedLymanAlpha()->intrinsicWavelengthRange());
    return range;
}

//////////////////////////////////////////////////////////////////////

double LyaSEDFamilyDecorator::specificLuminosity(double wavelength, const Array& parameters) const
{
    if (wavelength <= Constants::lambdaIon())
    {
        return (1. - conversionFraction()) * sedFamilyOriginal()->specificLuminosity(wavelength, parameters);
    }
    else
    {
        double luminosityLya = conversionFraction() * sedLymanAlpha()->specificLuminosity(wavelength);
        // calculating the ionizing luminosity is expensive, so we do this only when necessary
        if (luminosityLya > 0.)
        {
            Array lambdav, pv, Pv;  // the contents of these arrays is not used here
            luminosityLya *= sedFamilyOriginal()->cdf(lambdav, pv, Pv, _ionizingRange, parameters);
        }
        return sedFamilyOriginal()->specificLuminosity(wavelength, parameters) + luminosityLya;
    }
}

//////////////////////////////////////////////////////////////////////

double LyaSEDFamilyDecorator::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                                  const Array& parameters) const
{
    // determine the ionizing luminosity; the argument arrays lambdav, pv, Pv are used for temporary storage
    double ionizingLuminosity = sedFamilyOriginal()->cdf(lambdav, pv, Pv, _ionizingRange, parameters);

    // build a wavelength grid that includes wavelengths from both SEDs; the argument arrays are temporary storage
    vector<double> newlambdav;  // combined list of wavelengths
    sedFamilyOriginal()->cdf(lambdav, pv, Pv, wavelengthRange, parameters);
    for (double w : lambdav) newlambdav.push_back(w);
    sedLymanAlpha()->specificLuminosityArray(lambdav, pv, wavelengthRange);
    for (double w : lambdav) newlambdav.push_back(w);
    NR::unique(newlambdav);
    NR::assign(lambdav, newlambdav);

    // calculate the specific luminosity at each of these grid points
    size_t n = lambdav.size();
    pv.resize(n);
    for (size_t i = 0; i != n; ++i)
    {
        if (lambdav[i] <= Constants::lambdaIon())
            pv[i] = (1. - conversionFraction()) * sedFamilyOriginal()->specificLuminosity(lambdav[i], parameters);
        else
            pv[i] = sedFamilyOriginal()->specificLuminosity(lambdav[i], parameters)
                    + ionizingLuminosity * conversionFraction() * sedLymanAlpha()->specificLuminosity(lambdav[i]);
    }

    // calculate the cumulative distribution and normalization
    return NR::cdf2(true, lambdav, pv, Pv);
}

//////////////////////////////////////////////////////////////////////
