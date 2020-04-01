/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LyaUtils.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "Random.hpp"
#include "VoigtProfile.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    constexpr double c = Constants::c();              // speed of light in vacuum
    constexpr double kB = Constants::k();             // Boltzmann constant
    constexpr double mp = Constants::Mproton();       // proton mass
    constexpr double la = Constants::lambdaLya();     // central Lyman-alpha wavelength
    constexpr double Aa = Constants::EinsteinALya();  // Einstein A coefficient for Lyman-alpha transition
}

////////////////////////////////////////////////////////////////////

double LyaUtils::section(double lambda, double T)
{
    double vth = sqrt(2. * kB / mp * T);                 // thermal velocity for T
    double a = Aa * la / 4. / M_PI / vth;                // Voigt parameter
    double x = c / la * (la - lambda) / vth;             // dimensionless frequency
    double sigma0 = 3. * la * la * M_2_SQRTPI / 4. * a;  // cross section at line center
    return sigma0 * VoigtProfile::value(a, x);           // cross section at given x
}

////////////////////////////////////////////////////////////////////

std::pair<Vec, bool> LyaUtils::sampleAtomVelocity(double lambda, double T, double nH, double ds, Direction kin,
                                                  Configuration* config, Random* random)
{
    double vth = sqrt(2. * kB / mp * T);                 // thermal velocity for T
    double a = Aa * la / 4. / M_PI / vth;                // Voigt parameter
    double x = c / la * (la - lambda) / vth;             // dimensionless frequency
    double sigma0 = 3. * la * la * M_2_SQRTPI / 4. * a;  // cross section at line center

    // generate two directions that are orthogonal to each other and to the incoming photon packet direction
    Direction k1(1., 0., 0.);
    if (kin.kx() != 0. || kin.ky() != 0.)
    {
        k1 = Direction(kin.ky(), -kin.kx(), 0.);
        k1 /= k1.norm();
    }
    Direction k2(Vec::cross(k1, kin));

    // select the critical value of the dimensionless frequency depending on the acceleration scheme;
    // leaving the value at zero is equivalent to no acceleration
    double xcrit = 0.;
    switch (config->lyaAccelerationScheme())
    {
        case Configuration::LyaAccelerationScheme::None: break;
        case Configuration::LyaAccelerationScheme::Constant: xcrit = config->lyaAccelerationCriticalValue(); break;
        case Configuration::LyaAccelerationScheme::Laursen2009:
        {
            double atau = a * sigma0 * nH * ds;
            if (atau > 60)
                xcrit = 0.02 * std::exp(1.4 * pow(std::log(atau), 0.6));
            else if (atau > 1)
                xcrit = 0.02 * std::exp(0.6 * pow(std::log(atau), 1.2));
        }
        case Configuration::LyaAccelerationScheme::Smith2015:
        {
            double atau = a * sigma0 * nH * ds;
            if (atau >= 1) xcrit = 0.2 * std::cbrt(atau);
        }
    }

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

    // select the isotropic or the dipole phase function:
    // all wing events and 1/3 of core events are dipole, and the remaining 2/3 core events are isotropic,
    // where x=0.2 (in the atom frame) is used as cutoff between core and wing
    bool dipole = x > 0.2 || random->uniform() < 1. / 3.;

    // scale the atom velocity from dimensionless to regular units
    u *= vth;

    // return the atom velocity and the phase function choice
    return std::make_pair(u, dipole);
}

////////////////////////////////////////////////////////////////////

Range LyaUtils::relevantWavelengthRange(double vsmax, double vmmax, double nmax, double dmax)
{
    constexpr double tau = 1e-3;

    double vp_bulk = vsmax + vmmax;
    double vp_voigt = sqrt((3. * Aa * Aa * la * la * la * la) / (64. * M_PI * M_PI * M_PI * tau) * (nmax * dmax));
    double vp = vp_bulk + vp_voigt;

    return Range(la * (1 - vp / c), la * (1 + vp / c));
}

////////////////////////////////////////////////////////////////////
