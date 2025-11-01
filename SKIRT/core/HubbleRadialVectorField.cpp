/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "HubbleRadialVectorField.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

void HubbleRadialVectorField::setupSelfBefore()
{
    // verify property values
    if (_maxRadius < _turnoverRadius)
        throw FATALERROR("The maximum radius should not be smaller than the turnover radius");
}

////////////////////////////////////////////////////////////////////

int HubbleRadialVectorField::dimension() const
{
    return 3;
}

//////////////////////////////////////////////////////////////////////

Vec HubbleRadialVectorField::vector(Position bfr) const
{
    // get the radial distance
    double r = bfr.norm();

    // at the origin, or beyond rmax, immediately return a null vector
    if (r == 0 || r > maxRadius()) return Vec();

    // get a unit vector in the direction of the current position
    Vec u = bfr / r;

    // calculate the magnitude as a function of radial distance
    double v = 0.0;
    if (r <= turnoverRadius())
        v = r / turnoverRadius();
    else
        v = 1. - (r - turnoverRadius()) / (maxRadius() - turnoverRadius());

    // return a vector with proper magnitude and direction
    return v * u;
}

//////////////////////////////////////////////////////////////////////
