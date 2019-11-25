/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CylindricalVectorField.hpp"

//////////////////////////////////////////////////////////////////////

int CylindricalVectorField::dimension() const
{
    return 2;
}

//////////////////////////////////////////////////////////////////////

Vec CylindricalVectorField::vector(Position bfr) const
{
    Vec result(-bfr.y(), bfr.x(), 0.);
    double norm = result.norm();
    return norm > 0. ? result/norm : Vec();
}

//////////////////////////////////////////////////////////////////////
