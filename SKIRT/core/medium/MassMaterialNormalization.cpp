/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MassMaterialNormalization.hpp"
#include "MaterialMix.hpp"

////////////////////////////////////////////////////////////////////

std::pair<double, double> MassMaterialNormalization::numberAndMass(const Geometry* /*geom*/,
                                                                   const MaterialMix* mix) const
{
    return std::make_pair(_mass / mix->mass(), _mass);
}

////////////////////////////////////////////////////////////////////
