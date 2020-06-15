/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DipolePhaseFunction.hpp"
#include "Random.hpp"
#include "StokesVector.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // arbitrary constants determining the resolution for discretizing the scattering angles
    constexpr int Nphi = 361;  // phi from 0 to 2 pi, index f
    constexpr double dphi = 2 * M_PI / (Nphi - 1);
}

////////////////////////////////////////////////////////////////////

void DipolePhaseFunction::initialize(Random* random, bool includePolarization)
{
    // cache random number generator
    _random = random;

    if (includePolarization)
    {
        // create tables listing phi, phi/(2 pi), sin(2 phi) and 1-cos(2 phi) for a number of f indices
        _phiv.resize(Nphi);
        _phi1v.resize(Nphi);
        _phisv.resize(Nphi);
        _phicv.resize(Nphi);
        for (int f = 0; f < Nphi; f++)
        {
            double phi = f * dphi;
            _phiv[f] = phi;
            _phi1v[f] = phi / (2 * M_PI);
            _phisv[f] = sin(2 * phi);
            _phicv[f] = 1 - cos(2 * phi);
        }
    }
}

////////////////////////////////////////////////////////////////////

double DipolePhaseFunction::phaseFunctionValueForCosine(double costheta) const
{
    return 0.75 * (costheta * costheta + 1.);
}

////////////////////////////////////////////////////////////////////

double DipolePhaseFunction::generateCosineFromPhaseFunction() const
{
    double X = _random->uniform();
    double p = cbrt(4. * X - 2. + sqrt(16. * X * (X - 1.) + 5.));
    return p - 1. / p;
}

////////////////////////////////////////////////////////////////////

double DipolePhaseFunction::phaseFunctionValue(double theta, double phi, const StokesVector* sv) const
{
    double costheta = cos(theta);
    double S11 = costheta * costheta + 1.;
    double S12 = costheta * costheta - 1.;
    double polDegree = sv->linearPolarizationDegree();
    double polAngle = sv->polarizationAngle();
    return 0.75 * (S11 + polDegree * S12 * cos(2. * (phi - polAngle)));
}

////////////////////////////////////////////////////////////////////

std::pair<double, double> DipolePhaseFunction::generateAnglesFromPhaseFunction(const StokesVector* sv) const
{
    // get the theta cosine just as in the unpolarized case
    double costheta = generateCosineFromPhaseFunction();

    // construct and sample from the normalized cumulative distribution of phi for this wavelength and theta angle
    double polDegree = sv->linearPolarizationDegree();
    double polAngle = sv->polarizationAngle();
    double PF = polDegree * (costheta * costheta + 1.) / (costheta * costheta - 1.) / (4 * M_PI);
    double cos2polAngle = cos(2 * polAngle) * PF;
    double sin2polAngle = sin(2 * polAngle) * PF;
    double phi = _random->cdfLinLin(_phiv, _phi1v + cos2polAngle * _phisv + sin2polAngle * _phicv);

    // return the result
    return std::make_pair(acos(costheta), phi);
}

////////////////////////////////////////////////////////////////////

void DipolePhaseFunction::applyMueller(double theta, StokesVector* sv) const
{
    double costheta = cos(theta);
    sv->applyMueller(costheta * costheta + 1., costheta * costheta - 1., 2. * costheta, 0.);
}

////////////////////////////////////////////////////////////////////
