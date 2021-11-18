/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TabulatedWavelengthDistribution.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "SourceWavelengthRangeInterface.hpp"

//////////////////////////////////////////////////////////////////////

void TabulatedWavelengthDistribution::setupSelfBefore()
{
    WavelengthDistribution::setupSelfBefore();

    // obtain the wavelengths and probabilities
    Array inlambdav, inpv;
    getWavelengthsAndProbabilities(inlambdav, inpv);

    // reverse the arrays if needed to get the wavelengths in increasing order
    if (inlambdav.size() > 1 && inlambdav[0] > inlambdav[inlambdav.size() - 1])
    {
        std::reverse(begin(inlambdav), end(inlambdav));
        std::reverse(begin(inpv), end(inpv));
    }

    // determine the intersected wavelength range
    Range range(inlambdav[0], inlambdav[inlambdav.size() - 1]);
    range.intersect(interface<SourceWavelengthRangeInterface>()->wavelengthRange());
    if (range.empty()) throw FATALERROR("Wavelength distribution range does not overlap source wavelength range");

    // construct the regular and cumulative distributions in the intersected range
    NR::cdf<NR::interpolateLogLog>(_lambdav, _pv, _Pv, inlambdav, inpv, range);
}

//////////////////////////////////////////////////////////////////////

double TabulatedWavelengthDistribution::probability(double wavelength) const
{
    return NR::value<NR::interpolateLogLog>(wavelength, _lambdav, _pv);
}

//////////////////////////////////////////////////////////////////////

double TabulatedWavelengthDistribution::generateWavelength() const
{
    return random()->cdfLogLog(_lambdav, _pv, _Pv);
}

//////////////////////////////////////////////////////////////////////
