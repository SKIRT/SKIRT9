/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "HyperboloidGeometry.hpp"
#include "FatalError.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void HyperboloidGeometry::setupSelfBefore()
{
    AxGeometry::setupSelfBefore();

    // determine the maximal vertical extent of the hyperboloid
    _zmax = _D * cos(_Delta);

    // determine the radius of the top cross section
    _b = _D * sin(_Delta);

    // determine the imaginary axis of the hyperboloid
    _c = _a * _zmax / sqrt(_b * _b - _a * _a);

    // determine the normalization factor
    _A = 3. / 2. / M_PI / (_zmax * (2 * _a * _a + _b * _b) - 4. * pow(_rmin, 3));
}

//////////////////////////////////////////////////////////////////////

double HyperboloidGeometry::density(double R, double z) const
{
    double r = sqrt(R * R + z * z);
    double costheta = z / r;
    if (_rani)
    {
        double rminani = _rmin * sqrt(6. / 7. * fabs(costheta) * (2. * fabs(costheta) + 1));
        if (r <= rminani || r < _rcut) return 0.;
    }
    else
    {
        if (r <= _rmin) return 0.;
    }
    if (z > _zmax) return 0.;
    if (z < -_zmax) return 0.;
    if ((z >= 0.) && (z < _c / _a * sqrt(R * R - _a * _a))) return 0.;
    if ((z < 0.) && (z > -_c / _a * sqrt(R * R - _a * _a))) return 0.;
    return _A;
}

//////////////////////////////////////////////////////////////////////

Position HyperboloidGeometry::generatePosition() const
{
    while (true)
    {
        double R = _b * sqrt(random()->uniform());
        double z = _zmax * (2. * random()->uniform() - 1.);
        if (density(R, z))
        {
            double phi = 2. * M_PI * random()->uniform();
            return Position(R, phi, z, Position::CoordinateSystem::CYLINDRICAL);
        }
    }
}

//////////////////////////////////////////////////////////////////////

double HyperboloidGeometry::SigmaR() const
{
    return _A * (_a - _rmin);
}

//////////////////////////////////////////////////////////////////////

double HyperboloidGeometry::SigmaZ() const
{
    return 2. * _A * (_zmax - _rmin);
}

//////////////////////////////////////////////////////////////////////
