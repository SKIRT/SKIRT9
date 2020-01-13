/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "RangeGrainSizeDistribution.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

RangeGrainSizeDistribution::RangeGrainSizeDistribution(double minSize, double maxSize)
    : _minSize(minSize), _maxSize(maxSize)
{}

////////////////////////////////////////////////////////////////////

void RangeGrainSizeDistribution::setupSelfBefore()
{
    GrainSizeDistribution::setupSelfBefore();

    // verify the distribution range
    if (_maxSize <= _minSize) throw FATALERROR("Maximum grain size must be larger than minimum grain size");
}

////////////////////////////////////////////////////////////////////

double RangeGrainSizeDistribution::amin() const
{
    return _minSize;
}

////////////////////////////////////////////////////////////////////

double RangeGrainSizeDistribution::amax() const
{
    return _maxSize;
}

////////////////////////////////////////////////////////////////////
