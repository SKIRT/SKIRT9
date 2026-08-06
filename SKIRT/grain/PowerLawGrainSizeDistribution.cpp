/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PowerLawGrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

PowerLawGrainSizeDistribution::PowerLawGrainSizeDistribution(SimulationItem* parent, double minSize, double maxSize,
                                                             double exponent)
    : RangeGrainSizeDistribution(minSize, maxSize), _exponent(exponent)
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

double PowerLawGrainSizeDistribution::dnda(double a) const
{
    return pow(a, -_exponent);
}

////////////////////////////////////////////////////////////////////
