/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SepAxGeometry.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

Position SepAxGeometry::generatePosition() const
{
    // generate the random numbers in separate statements to guarantee evaluation order
    // (function arguments are evaluated in different order depending on the compiler)
    double R = randomCylRadius();
    double phi = 2.0 * M_PI * random()->uniform();
    double z = randomZ();
    return Position(R, phi, z, Position::CoordinateSystem::CYLINDRICAL);
}

//////////////////////////////////////////////////////////////////////
