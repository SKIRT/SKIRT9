/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StokesVector.hpp"

//////////////////////////////////////////////////////////////////////

StokesVector::StokesVector(double I, double Q, double U, double V, Direction n)
{
    setPolarized(I, Q, U, V, n);
}

//////////////////////////////////////////////////////////////////////

void StokesVector::setPolarized(double I, double Q, double U, double V, Direction n)
{
    if (I != 0.0 && !n.isNull())
    {
        _Q = Q / I;
        _U = U / I;
        _V = V / I;
        _normal = n;
        _polarized = true;
    }
    else
    {
        setUnpolarized();
    }
}

//////////////////////////////////////////////////////////////////////

void StokesVector::setPolarized(double I, double Q, double U, double V)
{
    if (I != 0.0)
    {
        _Q = Q / I;
        _U = U / I;
        _V = V / I;
        _polarized = true;
    }
    else
    {
        setUnpolarized();
    }
}

//////////////////////////////////////////////////////////////////////

void StokesVector::setPolarized(const StokesVector& polarization)
{
    *this = polarization;
}

//////////////////////////////////////////////////////////////////////

double StokesVector::totalPolarizationDegree() const
{
    return sqrt(_Q * _Q + _U * _U + _V * _V);
}

//////////////////////////////////////////////////////////////////////

double StokesVector::linearPolarizationDegree() const
{
    return sqrt(_Q * _Q + _U * _U);
}

//////////////////////////////////////////////////////////////////////

double StokesVector::polarizationAngle() const
{
    if (_U == 0 && _Q == 0)
        return 0;
    else
        return 0.5 * atan2(_U, _Q);
}

//////////////////////////////////////////////////////////////////////

void StokesVector::rotateStokes(double phi, Direction k)
{
    // if this is the first scattering: generate normal to scattering plane
    if (!_polarized)
    {
        double kx, ky, kz;
        k.cartesian(kx, ky, kz);
        // this is the Bianchi formula as used in Direction Random::direction(Direction bfk, double costheta)
        // with phi = 0 and theta = 90 deg.
        if (fabs(kz) > 0.99999)
        {
            _normal.set(1, 0, 0, false);
        }
        else
        {
            double nz = sqrt((1.0 - kz) * (1.0 + kz));
            double nx = -kx * kz / nz;
            double ny = -ky * kz / nz;
            _normal.set(nx, ny, nz, false);
        }
        _polarized = true;
    }
    else
    {
        // rotate the Q and U in the new reference frame
        double cos2phi = cos(2.0 * phi);
        double sin2phi = sin(2.0 * phi);
        double Q = cos2phi * _Q + sin2phi * _U;
        double U = -sin2phi * _Q + cos2phi * _U;
        _Q = Q;
        _U = U;
    }

    // rotate the stored scattering plane to obtain the new scattering plane using Rodrigues' rotation formula;
    // force (re)normalization to prevent degradation
    _normal = Direction(_normal * cos(phi) + Vec::cross(k, _normal) * sin(phi), true);
}

//////////////////////////////////////////////////////////////////////

double StokesVector::rotateIntoPlane(Direction k, Direction knew)
{
    // we want to rotate the reference into the plane, so we rotate the normal to be perpendicular to the plane.
    // generate the normal we want to rotate into
    Direction nNew(Vec::cross(k, knew), true);

    // if the given directions are parallel, simplify procedure
    if (nNew.isNull())
    {
        rotateStokes(0, k);
        return 0;
    }

    // if the StokesVector is unpolarized, simplify procedure
    if (!_polarized)
    {
        _normal = nNew;
        _polarized = true;
        return 0;
    }

    // calculate angle phi to rotate around
    double cosphi = Vec::dot(_normal, nNew);
    double sinphi = Vec::dot(Vec::cross(_normal, nNew), k);
    double phi = atan2(sinphi, cosphi);

    // rotate around the angle and return it
    rotateStokes(phi, k);
    return phi;
}

//////////////////////////////////////////////////////////////////////

void StokesVector::applyMueller(double S11, double S12, double S22, double S33, double S34, double S44)
{
    double I = S11 + S12 * _Q;
    double Q = S12 + S22 * _Q;
    double U = S33 * _U + S34 * _V;
    double V = -S34 * _U + S44 * _V;
    setPolarized(I, Q, U, V);
}

//////////////////////////////////////////////////////////////////////

void StokesVector::applyMueller(double S11, double S12, double S33, double S34)
{
    applyMueller(S11, S12, S11, S33, S34, S33);
}

//////////////////////////////////////////////////////////////////////
