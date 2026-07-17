/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NEBULARLINEEMISSION_HPP
#define NEBULARLINEEMISSION_HPP

#include "Basics.hpp"

//////////////////////////////////////////////////////////////////////

/** The NebularLineEmission class provides static data tables and functions for computing
    nebular emission from photoionized gas:

    - H recombination lines (Case B): Storey & Hummer (1995) P_B tables with Hui & Gnedin (1997)
      alpha_B. Lines: Lyman-alpha, Balmer (alpha through epsilon), Paschen (alpha, beta),
      Brackett alpha.
    - Metal forbidden lines: pre-tabulated collisional excitation rate coefficients q_col(T, n_e)
      from CHIANTI statistical equilibrium. Lines: [NII] 6548/6583, [OI] 6300/6364,
      [OII] 3727/3729, [OIII] 4363/4959/5007, [SII] 6716/6731.

    Implementation follows McClymont, Smith & Tacchella (2025).

    All internal component functions, coefficient tables, and helper constants are file-local to
    the implementation translation unit; only the public entry points are exposed. */
class NebularLineEmission final
{
public:
    // ============== Line definitions ==============

    // Line indices for the emission line array
    enum LineIndex {
        // H recombination lines (Storey & Hummer 1995 Case B)
        Lya = 0,   // Lyman-alpha   1215.67 A
        Ha,        // Balmer-alpha   6562.80 A
        Hb,        // Balmer-beta    4861.33 A
        Hg,        // Balmer-gamma   4340.46 A
        Hd,        // Balmer-delta   4101.73 A
        HeBalmer,  // Balmer-epsilon 3970.07 A
        Paa,       // Paschen-alpha  18751 A
        Pab,       // Paschen-beta   12818 A
        Bra,       // Brackett-alpha 40512 A

        // Metal forbidden lines (CHIANTI collisional excitation)
        NII6548,   // [NII] 6548
        NII6583,   // [NII] 6583
        OI6300,    // [OI] 6300
        OI6364,    // [OI] 6364
        OII3729,   // [OII] 3729  (2D 5/2)
        OII3726,   // [OII] 3726  (2D 3/2)
        OIII4363,  // [OIII] 4363
        OIII4959,  // [OIII] 4959
        OIII5007,  // [OIII] 5007
        SII6716,   // [SII] 6716
        SII6731,   // [SII] 6731
    };
    static constexpr int numLines = 20;

    // Rest-frame wavelengths [m]
    static constexpr double lineWavelengths[numLines] = {
        1215.670e-10,  // Lya
        6562.800e-10,  // Ha
        4861.330e-10,  // Hb
        4340.460e-10,  // Hg
        4101.730e-10,  // Hd
        3970.070e-10,  // He
        1.87510e-6,    // Paa
        1.28180e-6,    // Pab
        4.05120e-6,    // Bra
        6548.050e-10,  // [NII] 6548
        6583.460e-10,  // [NII] 6583
        6300.304e-10,  // [OI] 6300
        6363.776e-10,  // [OI] 6364
        3728.815e-10,  // [OII] 3729
        3726.032e-10,  // [OII] 3726
        4363.210e-10,  // [OIII] 4363
        4958.911e-10,  // [OIII] 4959
        5006.843e-10,  // [OIII] 5007
        6716.440e-10,  // [SII] 6716
        6730.810e-10,  // [SII] 6731
    };

    // Particle mass for Doppler broadening [kg]: proton mass for H, approximate for metals.
    static constexpr double lineMasses[numLines] = {
        1.67262192e-27,         // Lya (H)
        1.67262192e-27,         // Ha
        1.67262192e-27,         // Hb
        1.67262192e-27,         // Hg
        1.67262192e-27,         // Hd
        1.67262192e-27,         // He
        1.67262192e-27,         // Paa
        1.67262192e-27,         // Pab
        1.67262192e-27,         // Bra
        14.0 * 1.67262192e-27,  // NII
        14.0 * 1.67262192e-27,  // NII
        16.0 * 1.67262192e-27,  // OI
        16.0 * 1.67262192e-27,  // OI
        16.0 * 1.67262192e-27,  // OII
        16.0 * 1.67262192e-27,  // OII
        16.0 * 1.67262192e-27,  // OIII
        16.0 * 1.67262192e-27,  // OIII
        16.0 * 1.67262192e-27,  // OIII
        32.0 * 1.67262192e-27,  // SII
        32.0 * 1.67262192e-27,  // SII
    };

    // Ion stage index in PhotoIonizationSolver::ionFracs[] for carrier ion of each metal line
    // N: stages 12-19 (NI=12, NII=13, ...)
    // O: stages 20-28 (OI=20, OII=21, OIII=22, ...)
    // S: stages 62-66 (SI=62, SII=63, ...)
    static constexpr int lineCarrierIonIndex[numLines] = {
        -1, -1, -1, -1, -1, -1, -1, -1, -1,  // H lines: not from ion fractions
        13,                                  // [NII]  -> NII  = ion index 13
        13,                                  // [NII]  -> NII
        20,                                  // [OI]   -> OI   = ion index 20
        20,                                  // [OI]   -> OI
        21,                                  // [OII]  -> OII  = ion index 21
        21,                                  // [OII]  -> OII
        22,                                  // [OIII] -> OIII = ion index 22
        22,                                  // [OIII] -> OIII
        22,                                  // [OIII] -> OIII
        63,                                  // [SII]  -> SII  = ion index 63
        63,                                  // [SII]  -> SII
    };

    // Element index for each metal line (used to get abundance)
    // 0=C, 1=N, 2=O, 3=Ne, 4=Mg, 5=Si, 6=S, 7=Fe
    static constexpr int lineElementIndex[numLines] = {
        -1, -1, -1, -1, -1, -1, -1, -1, -1,  // H lines
        1,  1,                               // NII -> N
        2,  2,                               // OI -> O
        2,  2,                               // OII -> O
        2,  2,  2,                           // OIII -> O
        6,  6,                               // SII -> S
    };

    /** Compute H recombination line luminosity [W] in the photon-conserving form
        L = P_B(T, ne) * Gamma_HI * n_HI * V * h*nu_line. Assumes photoionization balance:
        each H-I photoionization eventually produces one Case B recombination, of which a
        fraction P_B emits in the requested line.
        lineIdx: LineIndex::Lya through LineIndex::Bra
        T: temperature [K]
        ne: electron density [cm^-3]
        gammaHI: H I photoionization rate per H atom [s^-1]
        nHI: neutral hydrogen number density [cm^-3]
        V_cm3: cell volume [cm^3]
        Returns luminosity in [W] (SI). */
    static double hydrogenLineLuminosity(int lineIdx, double T, double ne, double gammaHI, double nHI, double V_cm3);

    /** Compute metal forbidden line luminosity [W] for a given cell.
        lineIdx: LineIndex::NII6548 through LineIndex::SII6731
        T: temperature [K]
        ne: electron density [cm^-3]
        nIon: carrier ion number density [cm^-3]
        V_cm3: cell volume [cm^3]
        Returns luminosity in [W] (SI). */
    static double metalLineLuminosity(int lineIdx, double T, double ne, double nIon, double V_cm3);
};

//////////////////////////////////////////////////////////////////////

#endif
