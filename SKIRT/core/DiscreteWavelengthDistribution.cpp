/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DiscreteWavelengthDistribution.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "WavelengthRangeInterface.hpp"

//////////////////////////////////////////////////////////////////////

void DiscreteWavelengthDistribution::setupSelfBefore()
{
    WavelengthDistribution::setupSelfBefore();

    // get the source wavelength range
    Range range = interface<WavelengthRangeInterface>()->wavelengthRange();

    // obtain the wavelengths, pruning and sorting as needed
    vector<double> lambdav;
    for (double lambda : getWavelengths()) if (range.contains(lambda)) lambdav.push_back(lambda);
    NR::assign(_lambdav, lambdav);
    NR::sort(_lambdav);

    // verify that there is at least one wavelength and that the smallest one is positive
    size_t n = _lambdav.size();
    if (!n) throw FATALERROR("None of the wavelengths are in the source wavelength range");
    if (_lambdav[0] <= 0) throw FATALERROR("All wavelengths should be positive");

    // calculate the range borders
    _borderv.resize(2*n);
    double delta = _lambdav[0] * 1./1000.;
    for (size_t ell=0; ell!=n; ++ell)
    {
        _borderv[2*ell] = _lambdav[ell] - delta;
        _borderv[2*ell+1] = _lambdav[ell] + delta;
    }

    // verify that the ranges do not overlap
    if (!std::is_sorted(begin(_borderv), end(_borderv))) throw FATALERROR("Wavelengths are too close to each other");

    // calculate the probability within each range
    _probability = 1. / (2.*delta*n);
}

//////////////////////////////////////////////////////////////////////

double DiscreteWavelengthDistribution::probability(double wavelength) const
{
    // get the index of the phantom wavelength bin defined by the list of all K=N*2 borders
    //  0      => out of range on the left side
    //  K=N*2  => out of range on the right side
    //  even   => in between two ranges
    //  odd    => inside a range
    size_t index = std::upper_bound(begin(_borderv), end(_borderv), wavelength) - begin(_borderv);
    size_t inside = index & 1;    // lower bit
    return inside ? _probability : 0.;
}

//////////////////////////////////////////////////////////////////////

double DiscreteWavelengthDistribution::generateWavelength() const
{
    size_t index = static_cast<size_t>(random()->uniform() * _lambdav.size());
    return _lambdav[index];
}

//////////////////////////////////////////////////////////////////////
