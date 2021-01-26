/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "IsoThermalDiskGeometry.hpp"
#include "FatalError.hpp"
#include "Random.hpp"
#include "SpecialFunctions.hpp"

////////////////////////////////////////////////////////////////////

void IsoThermalDiskGeometry::setupSelfBefore()
{
    SepAxGeometry::setupSelfBefore();

    // calculate central density
    _rho0 = 1.0 / (4.0 * M_PI * _hz * _hR * _hR);
}

////////////////////////////////////////////////////////////////////

double IsoThermalDiskGeometry::density(double R, double z) const
{
    double sechz = 1.0 / cosh(z/_hz);
    return _rho0 * exp(-R / _hR) * sechz * sechz;
}

////////////////////////////////////////////////////////////////////

double IsoThermalDiskGeometry::randomCylRadius() const
{
    double X = random()->uniform();
    return _hR * (-1.0 - SpecialFunctions::LambertW1((X - 1.0) / M_E));
}

////////////////////////////////////////////////////////////////////

double IsoThermalDiskGeometry::randomZ() const
{
    double X = random()->uniform();
    return 0.5 * _hz * log(X / (1.0 - X));
}

//////////////////////////////////////////////////////////////////////

double IsoThermalDiskGeometry::SigmaR() const
{
    return _rho0 * _hR;
}

//////////////////////////////////////////////////////////////////////

double IsoThermalDiskGeometry::SigmaZ() const
{
    return 2.0 * _rho0 * _hz;
}

////////////////////////////////////////////////////////////////////
