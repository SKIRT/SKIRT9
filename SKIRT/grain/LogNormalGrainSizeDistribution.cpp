/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LogNormalGrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

LogNormalGrainSizeDistribution::LogNormalGrainSizeDistribution(SimulationItem* parent, double minSize, double maxSize,
                                                               double centroid, double width)
    : RangeGrainSizeDistribution(minSize, maxSize), _centroid(centroid), _width(width)
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

double LogNormalGrainSizeDistribution::dnda(double a) const
{
    double x = log(a / _centroid) / _width;
    return 1. / a * exp(-0.5 * x * x);
}

////////////////////////////////////////////////////////////////////
