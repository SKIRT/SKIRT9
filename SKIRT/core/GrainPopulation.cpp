/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GrainPopulation.hpp"

////////////////////////////////////////////////////////////////////

GrainPopulation::GrainPopulation(SimulationItem* parent,
                                 GrainComposition* composition, GrainSizeDistribution* sizeDistribution,
                                 int numSizes, NormalizationType normType, double normValue)
{
    parent->addChild(this);

    _composition = composition;
    _sizeDistribution = sizeDistribution;
    _numSizes = numSizes;

    switch (normType)
    {
    case NormalizationType::DustMassPerHydrogenAtom:
        _dustMassPerHydrogenAtom = normValue;
        break;
    case NormalizationType::DustMassPerHydrogenMass:
        _dustMassPerHydrogenMass = normValue;
        break;
    case NormalizationType::FactorOnSizeDistribution:
        _factorOnSizeDistribution = normValue;
        break;
    }

    setup();
}

////////////////////////////////////////////////////////////////////
