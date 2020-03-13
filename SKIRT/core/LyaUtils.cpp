/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LyaUtils.hpp"
#include "Constants.hpp"
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

double LyaUtils::sectionForDimlessFreq(double x, double T)
{
    double vth = sqrt(2. * kB * T / mp);               // thermal velocity for T
    double a = Aa * la / 4 / M_PI / vth;               // Voigt parameter
    double sigma0 = 3 * la * la * M_2_SQRTPI / 4 * a;  // cross section at line center
    return sigma0 * VoigtProfile::value(a, x);         // cross section at given x
}

////////////////////////////////////////////////////////////////////

double LyaUtils::sectionForWavelength(double lambda, double T)
{
    double vp = c * (lambda - la) / la;                // velocity shift for lambda
    double vth = sqrt(2. * kB * T / mp);               // thermal velocity for T
    double a = Aa * la / 4 / M_PI / vth;               // Voigt parameter
    double sigma0 = 3 * la * la * M_2_SQRTPI / 4 * a;  // cross section at line center
    return sigma0 * VoigtProfile::value(a, vp / vth);  // cross section at given x
}

////////////////////////////////////////////////////////////////////
