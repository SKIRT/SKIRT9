/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DonutGeometry.hpp"
#include "FatalError.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void DonutGeometry::setupSelfBefore()
{
    AxGeometry::setupSelfBefore();

    // verify property values
    if (_r0 > _R0)
        throw FATALERROR(
            "The radius in the equatorial plane should not be smaller than the radius in the azimuthal plane");

    // determine the normalization factor
    _A = 0.5 / M_PI / M_PI / _r0 / _r0 / _R0;

    // store the small circle radius squared
    _r02 = _r0 * _r0;
}

//////////////////////////////////////////////////////////////////////

double DonutGeometry::density(double R, double z) const
{
    // calculate the distance squared from the small circle centre
    double d2 = (_R0 - R) * (_R0 - R) + z * z;

    // the density is zero outside the ring torus geometry
    if (d2 >= _r02) return 0.0;

    // the density is uniform inside the ring torus geometry
    else
        return _A;
}

//////////////////////////////////////////////////////////////////////

Position DonutGeometry::generatePosition() const
{
    // generate a uniform azimuthal angle
    double phi = 2.0 * M_PI * random()->uniform();

    // generate a uniform position in the azimuthal plane circle
    double theta = 2.0 * M_PI * random()->uniform();
    double r = _r0 * sqrt(random()->uniform());

    // convert circle position to global cylindrical coordinates
    double z = r * sin(theta);
    double R = _R0 + r * cos(theta);

    // return position
    return Position(R, phi, z, Position::CoordinateSystem::CYLINDRICAL);
}

//////////////////////////////////////////////////////////////////////

double DonutGeometry::SigmaR() const
{
    return _A * 2 * _r0;
}

//////////////////////////////////////////////////////////////////////

double DonutGeometry::SigmaZ() const
{
    return 0.0;
}

//////////////////////////////////////////////////////////////////////
