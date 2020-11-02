/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ListGrainSizeDistribution.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

void ListGrainSizeDistribution::setupSelfBefore()
{
    // verify number of configured values
    if (_sizes.size() != _sizeDistributionValues.size())
        throw FATALERROR("Number of listed grain sizes does not match number of listed size distribution values");
    if (_sizes.size() < 2) throw FATALERROR("Number of listed grain sizes is less than 2");

    // copy the configured values into arrays (to allow the use of our regular interpolation function)
    _av = NR::array(_sizes);
    _dndav = NR::array(_sizeDistributionValues);
}

////////////////////////////////////////////////////////////////////

double ListGrainSizeDistribution::amin() const
{
    return *begin(_av);
}

////////////////////////////////////////////////////////////////////

double ListGrainSizeDistribution::amax() const
{
    return *(end(_av) - 1);
}

////////////////////////////////////////////////////////////////////

double ListGrainSizeDistribution::dnda(double a) const
{
    return NR::value<NR::interpolateLogLog>(a, _av, _dndav);
}

////////////////////////////////////////////////////////////////////
