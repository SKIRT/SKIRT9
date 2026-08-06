/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LyUtils.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "Random.hpp"
#include "VoigtProfile.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    constexpr double c = Constants::c();  // speed of light in vacuum
}

////////////////////////////////////////////////////////////////////

double LyUtils::section(double lambda, double center, double vth, double A, double a, double g)
{
    double x = (center - lambda) / lambda * c / vth;  // dimensionless frequency
    double sigma0 =
        g * center * center * center * M_2_SQRTPI / (16.0 * M_PI * vth) * A;  // cross section at line center
    return sigma0 * VoigtProfile::value(a, x);                                // cross section at given x
}

////////////////////////////////////////////////////////////////////

Vec LyUtils::sampleAtomVelocity(double lambda, double center, double vth, double a, double T, double nH, Direction kin,
                                Configuration* config, Random* random)
{
    double x = (center - lambda) / lambda * c / vth;  // dimensionless frequency

    // generate two directions that are orthogonal to each other and to the incoming photon packet direction
    Direction k1 = (kin.kx() || kin.ky()) ? Direction(kin.ky(), -kin.kx(), 0., true) : Direction(1., 0., 0., false);
    Direction k2(Vec::cross(k1, kin), false);

    // select the critical value of the dimensionless frequency depending on the acceleration scheme;
    // leaving the value at zero is equivalent to no acceleration
    double xcrit = 0.;
    switch (config->lyaAccelerationScheme())
    {
        case Configuration::LyaAccelerationScheme::None: break;
        case Configuration::LyaAccelerationScheme::Constant:
        {
            xcrit = config->lyaAccelerationStrength() * 3.;
            break;
        }
        case Configuration::LyaAccelerationScheme::Variable:
        {
            xcrit = config->lyaAccelerationStrength() * pow(nH / T, 1. / 6.);
            break;
        }
    }

    // apply the acceleration only to core scatterings, with xcrit defining the transition between core and wings
    if (abs(x) > xcrit) xcrit = 0.;

    // draw values for the components of the dimensionless atom velocity
    // parallel and orthogonal to the incoming photon packet
    double upar = VoigtProfile::sample(a, x, random);
    double radius = sqrt(xcrit * xcrit - std::log(random->uniform()));
    double angle = 2. * M_PI * random->uniform();
    double u1 = radius * cos(angle);
    double u2 = radius * sin(angle);

    // construct the dimensionless atom velocity from the direction vectors and the magnitudes
    Vec u = kin * upar + k1 * u1 + k2 * u2;

    // transform the dimensionless frequency into the rest frame of the atom
    x -= Vec::dot(u, kin);

    // scale the atom velocity from dimensionless to regular units
    u *= vth;

    // return the atom velocity and the phase function choice
    return u;
}

////////////////////////////////////////////////////////////////////

double LyUtils::shiftWavelength(double lambda, const Vec& vatom, const Direction& kin, const Direction& kout)
{
    return lambda / (1 - Vec::dot(kin, vatom) / c) * (1 - Vec::dot(kout, vatom) / c);
}

////////////////////////////////////////////////////////////////////
