/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "HirashitaLogNormalGrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

HirashitaLogNormalGrainSizeDistribution::HirashitaLogNormalGrainSizeDistribution(SimulationItem* parent, double minSize,
                                                                                 double maxSize, double centroid,
                                                                                 double width)
    : RangeGrainSizeDistribution(minSize, maxSize), _centroid(centroid), _width(width)
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

double HirashitaLogNormalGrainSizeDistribution::dnda(double a) const
{
    double x = log(a / _centroid) / _width;
    return 1. / pow(a, 4) * exp(-0.5 * x * x);
}

////////////////////////////////////////////////////////////////////
