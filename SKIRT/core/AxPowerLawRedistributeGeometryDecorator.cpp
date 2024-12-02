
/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AxPowerLawRedistributeGeometryDecorator.hpp"
#include "Log.hpp"

////////////////////////////////////////////////////////////////////

void AxPowerLawRedistributeGeometryDecorator::setupSelfAfter()
{
    RedistributeGeometryDecorator::setupSelfAfter();

    if (!_pR && !_pz) find<Log>()->warning("Both radial and vertical exponents of the axial power law are zero");
}

////////////////////////////////////////////////////////////////////

int AxPowerLawRedistributeGeometryDecorator::dimension() const
{
    return max(2, geometry()->dimension());
}

////////////////////////////////////////////////////////////////////

double AxPowerLawRedistributeGeometryDecorator::weight(Position bfr) const
{
    double R = bfr.cylRadius();
    double z = abs(bfr.z());
    if ((_pz && z < _z0) || (_pR && R < _R0)) return 0.;
    return pow(R, -_pR) * pow(z, -_pz);
}

////////////////////////////////////////////////////////////////////

double AxPowerLawRedistributeGeometryDecorator::maxWeight() const
{
    return weight(Position(_R0, 0., _z0, Position::CoordinateSystem::CYLINDRICAL));
}

////////////////////////////////////////////////////////////////////
