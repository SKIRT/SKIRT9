/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include "Basics.hpp"

//////////////////////////////////////////////////////////////////////

/** The Constants class offers static functions that return the value of physical constants or
    units used in astronomy. Other parts of the code should use these functions rather than
    defining the values locally. */
class Constants final
{
public:
    /** This function returns the speed of light \f$c\f$. */
    static double c();

    /** This function returns Planck's constant \f$h\f$. */
    static double h();

    /** This function returns the Boltzmann constant \f$k\f$. */
    static double k();

    /** This function returns the Avogadro constant \f$N_\text{A}\f$, i.e. the number of atoms or
        molecules in one mole of a given substance. */
    static double NA();

    /** This function returns the distance of one astronomical unit. */
    static double AU();

    /** This function returns the distance of one parsec. */
    static double pc();

    /** This function returns the rest mass of the proton. */
    static double Mproton();

    /** This function returns the rest mass of the electron. */
    static double Melectron();

    /** This function returns the mass of the hydrogen atom. */
    static double MHatom();

    /** This function returns the solar mass. */
    static double Msun();

    /** This function returns the solar bolometric luminosity without neutrino radiation. */
    static double Lsun();

    /** This function returns the average temperature of the cosmic microwave background (CMB) in
        the Local Universe. */
    static double Tcmb();

    /** This function returns the ionization wavelength of atomic hydrogen. */
    static double lambdaIon();

    /** This function returns the central wavelength of the hydrogen Lyman-alpha transition. */
    static double lambdaLya();

    /** This function returns the wavelength of the center of the optical \f$V\f$ band. */
    static double lambdaV();

    /** This function returns the total (i.e. absorption and scattering) dust opacity
        \f$\kappa_{\text{V}}\f$ in the V-band for a Draine & Li (2007) dust mix with
        characteristics corresponding to the InterstellarDistMix class. */
    static double kappaV();

    /** This function returns the Thomson cross section for an electron. */
    static double sigmaThomson();

    /** This function returns the length of a Julian year in seconds. A Julian year has
        exactly 365.25 days of exactly 86400 seconds each. */
    static double year();
};

//////////////////////////////////////////////////////////////////////

#endif
