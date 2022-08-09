/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "NumberMaterialNormalization.hpp"
#include "MaterialMix.hpp"

////////////////////////////////////////////////////////////////////

std::pair<double, double> NumberMaterialNormalization::numberAndMass(const Geometry* /*geom*/,
                                                                     const MaterialMix* mix) const
{
    return std::make_pair(_number, _number * mix->mass());
}

////////////////////////////////////////////////////////////////////
