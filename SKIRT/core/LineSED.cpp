/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LineSED.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void LineSED::setupSelfBefore()
{
    SED::setupSelfBefore();

    // obtain the line wavelengths and luminosities from the subclass
    getWavelengthsAndLuminosities(_inlambdav, _inLv);

    // verify that there is at least one intrinsic line
    size_t numLines = _inlambdav.size();
    if (!numLines) throw FATALERROR("Line SED must have at least one wavelength/luminosity pair");

    // copy the lines inside the source wavelength range to temporary vectors, and then to data member arrays
    auto range = find<Configuration>()->sourceWavelengthRange();
    vector<double> lambdav;
    vector<double> Lv;
    for (size_t i = 0; i != numLines; ++i)
    {
        if (range.containsFuzzy(_inlambdav[i]))
        {
            lambdav.push_back(_inlambdav[i]);
            Lv.push_back(_inLv[i]);
        }
    }
    NR::assign(_lambdav, lambdav);
    NR::assign(_Lv, Lv);

    // verify that there is at least one line in the source wavelength range
    if (!_lambdav.size()) throw FATALERROR("Line SED must have at least one line in the source wavelength range");

    // construct the cumulative distribution and normalize the luminosities
    double norm = NR::cdf(_Pv, _Lv);
    _inLv /= norm;
    _Lv /= norm;
}

////////////////////////////////////////////////////////////////////

Range LineSED::intrinsicWavelengthRange() const
{
    // avoid empty range if there is just a single line
    return Range(_inlambdav.min() * (1 - 1e-6), _inlambdav.max() * (1 + 1e-6));
}

////////////////////////////////////////////////////////////////////

double LineSED::integratedLuminosity(const Range& wavelengthRange) const
{
    double sum = 0.;
    size_t numLines = _inlambdav.size();
    for (size_t i = 0; i != numLines; ++i)
    {
        if (wavelengthRange.containsFuzzy(_inlambdav[i])) sum += _inLv[i];
    }
    return sum;
}

////////////////////////////////////////////////////////////////////

double LineSED::generateWavelength() const
{
    // select a line with a probability corresponding to its luminosity contribution
    int i = NR::locateClip(_Pv, random()->uniform());
    return _lambdav[i];
}

////////////////////////////////////////////////////////////////////
