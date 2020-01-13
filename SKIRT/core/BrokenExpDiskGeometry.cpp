/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BrokenExpDiskGeometry.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void BrokenExpDiskGeometry::setupSelfBefore()
{
    SepAxGeometry::setupSelfBefore();

    _beta = 1.0 / _s * (_hout / _hinn - 1.0);

    // create a radial array with the cumulative mass distribution. Use Ninn points in the
    // range between 0 and _Rb, and Nout points between _Rb and the outermost radius,
    // which we choose to be _Rb + 10*_hout.

    int Ninn = 200;
    int Nout = 400;
    int N = Ninn + Nout;
    _Rv.resize(N + 1);
    _Xv.resize(N + 1);
    double dRinn = _Rb / Ninn;
    for (int i = 0; i < Ninn; i++) _Rv[i] = i * dRinn;
    double dRout = 10 * _hout / Nout;
    for (int i = Ninn; i <= N; i++) _Rv[i] = _Rb + (i - Ninn) * dRout;
    _Xv[0] = 0.0;
    double intprev = 0.0;
    for (int i = 1; i <= N; i++)
    {
        double RL = _Rv[i - 1];
        double RR = _Rv[i];
        double intL = intprev;
        double intR = radialDensity(RR) * RR;
        intprev = intR;
        _Xv[i] = _Xv[i - 1] + 0.5 * (RR - RL) * (intL + intR);
    }
    double IR = _Xv[N - 1];
    _Xv /= IR;

    // calculate _rho0;

    _rho0 = 1.0 / 4.0 / M_PI / _hz / IR;

    // calculate the radial surface density

    intprev = 0.0;
    double sum = 0.0;
    for (int i = 1; i <= N; i++)
    {
        double RL = _Rv[i - 1];
        double RR = _Rv[i];
        double intL = intprev;
        double intR = radialDensity(RR);
        intprev = intR;
        sum += 0.5 * (RR - RL) * (intL + intR);
    }
    _SigmaR = sum * _rho0;
}

////////////////////////////////////////////////////////////////////

double BrokenExpDiskGeometry::radialDensity(double R) const
{
    return exp(-R / _hinn) * pow(1.0 + exp(_s * (R - _Rb) / _hout), _beta);
}

////////////////////////////////////////////////////////////////////

double BrokenExpDiskGeometry::density(double R, double z) const
{
    return _rho0 * exp(-abs(z) / _hz) * radialDensity(R);
}

////////////////////////////////////////////////////////////////////

double BrokenExpDiskGeometry::randomCylRadius() const
{
    return random()->cdfLinLin(_Rv, _Xv);
}

////////////////////////////////////////////////////////////////////

double BrokenExpDiskGeometry::randomZ() const
{
    double z, XX;
    XX = random()->uniform();
    z = (XX <= 0.5) ? _hz * log(2.0 * XX) : -_hz * log(2.0 * (1.0 - XX));
    return z;
}

//////////////////////////////////////////////////////////////////////

double BrokenExpDiskGeometry::SigmaR() const
{
    return _SigmaR;
}

//////////////////////////////////////////////////////////////////////

double BrokenExpDiskGeometry::SigmaZ() const
{
    return 2.0 * _rho0 * _hz * radialDensity(0.0);
}

//////////////////////////////////////////////////////////////////////
