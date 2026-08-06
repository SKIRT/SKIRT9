/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpiralStructureGeometryDecorator.hpp"
#include "Random.hpp"
#include "SpecialFunctions.hpp"

//////////////////////////////////////////////////////////////////////

void SpiralStructureGeometryDecorator::setupSelfBefore()
{
    GenGeometry::setupSelfBefore();

    // cache frequently used values
    _tanp = tan(_p);
    _cn = sqrt(M_PI) * SpecialFunctions::gamma(_N + 1.0) / SpecialFunctions::gamma(_N + 0.5);
    _c = 1.0 + (_cn - 1.0) * _w;
}

////////////////////////////////////////////////////////////////////

double SpiralStructureGeometryDecorator::density(Position bfr) const
{
    double R, phi, z;
    bfr.cylindrical(R, phi, z);
    return _geometry->density(R, z) * perturbation(R, phi);
}

////////////////////////////////////////////////////////////////////

Position SpiralStructureGeometryDecorator::generatePosition() const
{
    Position bfr = _geometry->generatePosition();
    double R, dummyphi, z;
    bfr.cylindrical(R, dummyphi, z);
    double c = 1.0 + (_cn - 1.0) * _w;
    double phi, t;
    do
    {
        phi = 2.0 * M_PI * random()->uniform();
        t = random()->uniform() * c / perturbation(R, phi);
    } while (t > 1);
    return Position(R, phi, z, Position::CoordinateSystem::CYLINDRICAL);
}

////////////////////////////////////////////////////////////////////

double SpiralStructureGeometryDecorator::SigmaX() const
{
    return _geometry->SigmaX();
}

////////////////////////////////////////////////////////////////////

double SpiralStructureGeometryDecorator::SigmaY() const
{
    return _geometry->SigmaY();
}

////////////////////////////////////////////////////////////////////

double SpiralStructureGeometryDecorator::SigmaZ() const
{
    return _geometry->SigmaZ();
}

////////////////////////////////////////////////////////////////////

double SpiralStructureGeometryDecorator::perturbation(double R, double phi) const
{
    double gamma = log(R / _R0) / _tanp + _phi0 + 0.5 * M_PI / _m;
    return (1.0 - _w) + _w * _cn * pow(sin(0.5 * _m * (gamma - phi)), 2 * _N);
}

////////////////////////////////////////////////////////////////////
