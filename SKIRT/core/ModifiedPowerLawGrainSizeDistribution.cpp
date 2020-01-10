/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ModifiedPowerLawGrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

ModifiedPowerLawGrainSizeDistribution::ModifiedPowerLawGrainSizeDistribution(
    SimulationItem* parent, double minSize, double maxSize, double powerLawIndex, double turnOffPoint,
    double scaleExponentialDecay, double exponentExponentialDecay, double scaleCurvature, double strengthCurvature,
    double exponentCurvature)

    : RangeGrainSizeDistribution(minSize, maxSize), _powerLawIndex(powerLawIndex), _turnOffPoint{turnOffPoint},
      _scaleExponentialDecay{scaleExponentialDecay}, _exponentExponentialDecay{exponentExponentialDecay},
      _scaleCurvature{scaleCurvature}, _strengthCurvature{strengthCurvature}, _exponentCurvature{exponentCurvature}
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

double ModifiedPowerLawGrainSizeDistribution::dnda(double a) const
{
    return pow(a, _alpha) * pow(1. + fabs(_zeta) * pow(a / _au, _eta), _zeta >= 0 ? 1 : -1)  // curvature
           * (a <= _at ? 1. : exp(-pow((a - _at) / _ac, _gamma)));                           // exponential decay
}

////////////////////////////////////////////////////////////////////
