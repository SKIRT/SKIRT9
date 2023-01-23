/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AnnulusGeometry.hpp"
#include "FatalError.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void AnnulusGeometry::setupSelfBefore()
{
    SepAxGeometry::setupSelfBefore();

    // verify property values
    if (_Rmin >= _Rmax) throw FATALERROR("The radius of the central cavity should be smaller than the outer radius");

    // calculate the uniform density
    _rho0 = 1.0 / ((_Rmax * _Rmax - _Rmin * _Rmin) * M_PI * _h);
}

////////////////////////////////////////////////////////////////////

double AnnulusGeometry::density(double R, double z) const
{
    double absz = fabs(z);
    if (R > _Rmax)
        return 0.0;
    else if (absz > _h / 2.0)
        return 0.0;
    else if (R < _Rmin)
        return 0.0;
    return _rho0;
}

////////////////////////////////////////////////////////////////////

double AnnulusGeometry::randomCylRadius() const
{
    return sqrt(_Rmin * _Rmin + (_Rmax - _Rmin) * (_Rmax + _Rmin) * random()->uniform());
}

////////////////////////////////////////////////////////////////////

double AnnulusGeometry::randomZ() const
{
    double X = random()->uniform();
    return (X - 0.5) * _h;
}

//////////////////////////////////////////////////////////////////////

double AnnulusGeometry::SigmaR() const
{
    return _rho0 * (_Rmax - _Rmin);
}

//////////////////////////////////////////////////////////////////////

double AnnulusGeometry::SigmaZ() const
{
    if (_Rmin > 0.0)
        return 0.0;
    else
        return _rho0 * _h;
}

////////////////////////////////////////////////////////////////////
