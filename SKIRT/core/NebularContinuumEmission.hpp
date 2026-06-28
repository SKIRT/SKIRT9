/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NEBULARCONTINUUMEMISSION_HPP
#define NEBULARCONTINUUMEMISSION_HPP

#include "Basics.hpp"

//////////////////////////////////////////////////////////////////////

/** The NebularContinuumEmission namespace provides functions for computing the nebular continuum
    emission spectrum (free-bound, free-free, and two-photon) from ionized hydrogen and helium.

    The implementation follows McClymont, Smith & Tacchella (2025), which is based on
    NEBULAR (Schirmer 2016):
    - Free-bound: Ercolano & Storey (2006) recombination coefficients for HI, HeI, HeII
    - Free-free: van Hoof et al. (2014) thermally averaged Gaunt factors with 3rd-order
      Lagrange interpolation on an 81×146 grid in log(γ²) × log(u)
    - Two-photon: Nussbaumer & Schmutz (1984) spectral shape with Hummer & Storey (1987)
      effective recombination rates and Pengelly & Seaton (1964) collisional de-excitation

    All emissivity functions return values in CGS units [erg s⁻¹ Hz⁻¹ cm⁻³] and must be
    multiplied by ne × nion × V to obtain the total luminosity per wavelength bin.

    The full continuum emissivity per unit wavelength [W/m] for a cell is:
        j_lambda = (j_ff + j_fb_HI + j_fb_HeI + j_fb_HeII + j_2p_HI + j_2p_HeII)
                    × (c / lambda²) × ne × nion × V × 1e-7
    where the (c/lambda²) converts from per-Hz to per-wavelength and 1e-7 converts erg/s to W.
*/
namespace NebularContinuumEmission
{
    // ===== Physical constants (CGS) =====
    constexpr double h_cgs = 6.62607015e-27;                      // Planck constant [erg s]
    constexpr double c_cgs = 2.99792458e10;                       // Speed of light [cm/s]
    constexpr double kB_cgs = 1.380649e-16;                       // Boltzmann constant [erg/K]
    constexpr double eV_cgs = 1.602176634e-12;                    // eV in erg
    constexpr double Ryd_cgs = 13.605693122994 * eV_cgs;          // Rydberg energy [erg]
    constexpr double Ryd_K = Ryd_cgs / kB_cgs;                    // Rydberg energy in Kelvin
    constexpr double angstrom = 1e-8;                             // Angstrom in cm
    constexpr double E_w = h_cgs * c_cgs / (angstrom * Ryd_cgs);  // hc/(Å Ry) [Å Ry / Ry = dimensionless]
    constexpr double mHCgs = 1.6735575e-24;                       // Hydrogen mass [g]
    constexpr double mHe_cgs = 6.6464764e-24;                     // Helium-4 mass [g]
    constexpr double meCgs = 9.1093837e-28;                       // Electron mass [g]
    constexpr double ee_cgs = 4.8032047e-10;                      // Elementary charge [esu]

    // Two-photon transition frequencies and rates
    constexpr double nu_12_HI = c_cgs * 82258.9544;                // HI 2s→1s 2p transition frequency [Hz]
    constexpr double nu_12_HeII = c_cgs * 329179.7623;             // HeII 2s→1s 2p frequency [Hz]
    constexpr double Z6_HeII = 64.0 * 1.097373156e5 / 1.096788e5;  // Z^6 (R_Z/R_H) for HeII
    constexpr double A_2q_HI = 8.2249;                             // HI 2s→1s two-photon rate [s⁻¹]
    constexpr double A_2q_HeII = A_2q_HI * Z6_HeII;                // HeII 2s→1s two-photon rate [s⁻¹]

    // ===== Free-bound continuum =====

    /** Compute the HI free-bound emission coefficient at wavelength w [Angstrom] and temperature T [K].
        Returns gamma_fb [dimensionless] such that j_nu = gamma_fb * ne * nHII [erg/s/Hz/cm³].
        Uses Ercolano & Storey (2006) tables with bilinear interpolation. */
    double HI_fb(double T, double w_angstrom);

    /** Compute the HeI free-bound emission coefficient (same convention as HI_fb). */
    double HeI_fb(double T, double w_angstrom);

    /** Compute the HeII free-bound emission coefficient (same convention as HI_fb). */
    double HeII_fb(double T, double w_angstrom);

    // ===== Free-free continuum =====

    /** Evaluate the thermally averaged free-free Gaunt factor g_ff(gamma^2, u)
        where gamma^2 = Z^2 Ry / (kB T) and u = h nu / (kB T).
        Uses van Hoof et al. (2014) 81×146 table with 3rd-order Lagrange interpolation.
        Arguments are log10(gamma^2) and log10(u). Returns 0 if outside grid. */
    double gauntFactor(double log_gamma2, double log_u);

    /** Compute the free-free emission coefficient at wavelength w [Angstrom] for charge Z and temperature T [K].
        Returns j_ff [erg s⁻¹ Hz⁻¹ cm³] (per ne per nion).
        j_ff = ff_prefactor * Z² * g_ff * exp(-u) / sqrt(T) / w²
        where u = hc / (lambda kB T). */
    double freeFree(double T, double w_angstrom, int Z);

    // ===== Two-photon continuum =====

    /** Compute the HI two-photon emission rate alpha_eff * A_2q / (A_2q + q_coll)
        at temperature T [K] and electron density ne [cm⁻³].
        Must be multiplied by g_nu(w) to get the frequency-dependent emissivity.
        Returns rate [cm³ s⁻¹] (per ne per nHII, then multiply by g_nu). */
    double HI_2p_rate(double T, double ne, double nHII);

    /** Compute the HeII two-photon emission rate (same convention as HI_2p_rate). */
    double HeII_2p_rate(double T, double ne, double nHeIII);

    /** Compute the HI two-photon spectral shape g_nu(w) at wavelength w [Angstrom].
        g_nu = h * (nu/nu_12) * A_nu / A_2q following Nussbaumer & Schmutz (1984).
        Units: [erg]. Total two-photon emissivity = rate(T,ne,nHII) * g_nu(w) * ne * nHII. */
    double HI_2p_gnu(double w_angstrom);

    /** Compute the HeII two-photon spectral shape g_nu(w) (same convention). */
    double HeII_2p_gnu(double w_angstrom);

    // ===== Full continuum spectrum =====

    /** Compute the full nebular continuum emissivity at wavelength lambda [m] for a cell with:
        T [K], ne [cm⁻³], nHII [cm⁻³], nHeII [cm⁻³], nHeIII [cm⁻³], V [cm³].
        Returns luminosity per unit wavelength [W/m] at this wavelength.
        This is what SKIRT's emissionSpectrum() should return for each wavelength bin. */
    double continuumLuminosity(double lambda_m, double T, double ne, double nHII, double nHeII, double nHeIII,
                               double V_cm3);

}  // namespace NebularContinuumEmission

//////////////////////////////////////////////////////////////////////

#endif
