/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TabulatedWavelengthDistribution.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "WavelengthRangeInterface.hpp"

//////////////////////////////////////////////////////////////////////

void TabulatedWavelengthDistribution::setupSelfBefore()
{
    WavelengthDistribution::setupSelfBefore();

    // obtain the wavelengths and probabilities
    Array inlambdav, inpv;
    getWavelengthsAndProbabilities(inlambdav, inpv);

    // determine the intersected wavelength range
    Range range(inlambdav[0], inlambdav[inlambdav.size()-1]);
    range.intersect(interface<WavelengthRangeInterface>()->wavelengthRange());
    if (range.empty()) throw FATALERROR("Wavelength distribution range does not overlap source wavelength range");

    // construct the regular and cumulative distributions in the intersected range
    NR::cdf<NR::interpolateLogLog>(_lambdav, _pv, _Pv, inlambdav, inpv, range);
}

//////////////////////////////////////////////////////////////////////

double TabulatedWavelengthDistribution::probability(double wavelength) const
{
    int i = NR::locateFail(_lambdav, wavelength);
    if (i < 0) return 0.;
    return NR::interpolateLogLog(wavelength, _lambdav[i], _lambdav[i+1], _pv[i], _pv[i+1]);
}

//////////////////////////////////////////////////////////////////////

double TabulatedWavelengthDistribution::generateWavelength() const
{
    return random()->cdfLogLog(_lambdav, _pv, _Pv);
}

//////////////////////////////////////////////////////////////////////
