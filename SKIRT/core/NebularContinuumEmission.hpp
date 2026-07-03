/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NEBULARCONTINUUMEMISSION_HPP
#define NEBULARCONTINUUMEMISSION_HPP

#include "Basics.hpp"

//////////////////////////////////////////////////////////////////////

/** The NebularContinuumEmission class provides functions for computing the nebular continuum
    emission spectrum (free-bound, free-free, and two-photon) from ionized hydrogen and helium.

    The implementation follows McClymont, Smith & Tacchella (2025), which is based on
    NEBULAR (Schirmer 2016):
    - Free-bound: Ercolano & Storey (2006) recombination coefficients for HI, HeI, HeII
    - Free-free: van Hoof et al. (2014) thermally averaged Gaunt factors with 3rd-order
      Lagrange interpolation on an 81x146 grid in log(gamma^2) x log(u)
    - Two-photon: Nussbaumer & Schmutz (1984) spectral shape with Hummer & Storey (1987)
      effective recombination rates and Pengelly & Seaton (1964) collisional de-excitation

    All internal component functions and coefficient tables are file-local to the
    implementation translation unit; only the full-spectrum entry point is exposed. */
class NebularContinuumEmission final
{
public:
    /** Compute the full nebular continuum emissivity at wavelength lambda [m] for a cell with:
        T [K], ne [cm^-3], nHII [cm^-3], nHeII [cm^-3], nHeIII [cm^-3], V [cm^3].
        Returns luminosity per unit wavelength [W/m] at this wavelength.
        This is what SKIRT's emissionSpectrum() should return for each wavelength bin. */
    static double continuumLuminosity(double lambda_m, double T, double ne, double nHII, double nHeII, double nHeIII,
                                      double V_cm3);
};

//////////////////////////////////////////////////////////////////////

#endif
