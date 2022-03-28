/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include "Basics.hpp"

//////////////////////////////////////////////////////////////////////

/** The Constants namespace offers functions that return the value of physical constants or units
    used in astronomy. Other parts of the code should use these functions rather than defining the
    values locally. All functions are defined inline as \c constexpr to avoid runtime overhead. */
namespace Constants
{
    /** This function returns the speed of light \f$c\f$. */
    constexpr double c() { return 2.99792458e8; }

    /** This function returns Planck's constant \f$h\f$. */
    constexpr double h() { return 6.62606957e-34; }

    /** This function returns the Boltzmann constant \f$k\f$. */
    constexpr double k() { return 1.3806488e-23; }

    /** This function returns the Avogadro constant \f$N_\text{A}\f$, i.e. the number of atoms or
        molecules in one mole of a given substance. */
    constexpr double NA() { return 6.02214129e23; }

    /** This function returns the distance of one astronomical unit. */
    constexpr double AU() { return 1.49597871e11; }

    /** This function returns the distance of one parsec. */
    constexpr double pc() { return 3.08567758e16; }

    /** This function returns the mass of an atomic mass unit (1/12 of the mass of an unbound
        neutral atom of carbon-12 in its nuclear and electronic ground state and at rest). */
    constexpr double amu() { return 1.6605390666e-27; }

    /** This function returns the rest mass of the proton. */
    constexpr double Mproton() { return 1.67262178e-27; }

    /** This function returns the rest mass of the electron. */
    constexpr double Melectron() { return 9.10938215e-31; }

    /** This function returns the electric charge of an electron in coulombs, which is equivalent
        to the value of 1 eV in joules. */
    constexpr double Qelectron() { return 1.602176634e-19; }

    /** This function returns the solar mass. */
    constexpr double Msun() { return 1.9891e30; }

    /** This function returns the solar bolometric luminosity without neutrino radiation. */
    constexpr double Lsun() { return 3.839e26; }

    /** This function returns the average temperature of the cosmic microwave background (CMB) in
        the Local Universe. */
    constexpr double Tcmb() { return 2.725; }

    /** This function returns the ionization wavelength of atomic hydrogen. */
    constexpr double lambdaIon() { return 911.75e-10; }

    /** This function returns the central wavelength of the hydrogen Lyman-alpha transition. */
    constexpr double lambdaLya() { return 1215.67e-10; }

    /** This function returns the Einstein A-coefficient of the Lyman-alpha transition. */
    constexpr double EinsteinALya() { return 6.25e8; }

    /** This function returns the Thomson cross section for an electron. */
    constexpr double sigmaThomson() { return 6.6524587158e-29; }

    /** This function returns the length of a Julian year in seconds. A Julian year has
        exactly 365.25 days of exactly 86400 seconds each. */
    constexpr double year() { return 31557600.; }
}

//////////////////////////////////////////////////////////////////////

#endif
