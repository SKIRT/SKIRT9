/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef VERNERCROSSSECTIONS_HPP
#define VERNERCROSSSECTIONS_HPP

#include "Basics.hpp"

//////////////////////////////////////////////////////////////////////

/** The VernerCrossSections class provides photoionization cross-sections for all relevant
    ions of H, He, C, N, O, Ne, Mg, Si, S, and Fe, using the analytic fitting formulae of
    Verner et al. (1996, ApJ, 465, 487). Each function returns the cross-section in cm^2 for a
    given photon energy in eV, or zero if the energy is below the ionization threshold or above
    the fitted range.

    The individual per-ion fitting formulae are file-local to the implementation translation
    unit; only the dispatch entry points, ionization potentials, and index-range constants are
    exposed here. */
class VernerCrossSections final
{
public:
    // ---- ionization potentials [eV] ----

    static constexpr double HI_eV = 13.6;
    static constexpr double HeI_eV = 24.59;
    static constexpr double HeII_eV = 54.42;

    static constexpr double CI_eV = 11.26;
    static constexpr double CII_eV = 24.38;
    static constexpr double CIII_eV = 47.89;
    static constexpr double CIV_eV = 64.49;
    static constexpr double CV_eV = 392.1;
    static constexpr double CVI_eV = 490.;

    static constexpr double NI_eV = 14.53;
    static constexpr double NII_eV = 29.6;
    static constexpr double NIII_eV = 47.45;
    static constexpr double NIV_eV = 77.47;
    static constexpr double NV_eV = 97.89;
    static constexpr double NVI_eV = 552.1;
    static constexpr double NVII_eV = 667.1;

    static constexpr double OI_eV = 13.62;
    static constexpr double OII_eV = 35.12;
    static constexpr double OIII_eV = 54.94;
    static constexpr double OIV_eV = 77.41;
    static constexpr double OV_eV = 113.9;
    static constexpr double OVI_eV = 138.1;
    static constexpr double OVII_eV = 739.2;
    static constexpr double OVIII_eV = 871.4;

    static constexpr double NeI_eV = 21.56;
    static constexpr double NeII_eV = 40.96;
    static constexpr double NeIII_eV = 63.46;
    static constexpr double NeIV_eV = 97.12;
    static constexpr double NeV_eV = 126.2;
    static constexpr double NeVI_eV = 157.9;
    static constexpr double NeVII_eV = 207.3;
    static constexpr double NeVIII_eV = 239.1;

    static constexpr double MgI_eV = 7.646;
    static constexpr double MgII_eV = 15.04;
    static constexpr double MgIII_eV = 80.14;
    static constexpr double MgIV_eV = 109.3;
    static constexpr double MgV_eV = 141.3;
    static constexpr double MgVI_eV = 186.5;
    static constexpr double MgVII_eV = 224.9;
    static constexpr double MgVIII_eV = 266.;
    static constexpr double MgIX_eV = 328.2;
    static constexpr double MgX_eV = 367.5;

    static constexpr double SiI_eV = 8.152;
    static constexpr double SiII_eV = 16.35;
    static constexpr double SiIII_eV = 33.49;
    static constexpr double SiIV_eV = 45.14;
    static constexpr double SiV_eV = 166.8;
    static constexpr double SiVI_eV = 205.1;
    static constexpr double SiVII_eV = 246.5;
    static constexpr double SiVIII_eV = 303.2;
    static constexpr double SiIX_eV = 351.1;
    static constexpr double SiX_eV = 401.4;
    static constexpr double SiXI_eV = 476.1;
    static constexpr double SiXII_eV = 523.5;

    static constexpr double SI_eV = 10.36;
    static constexpr double SII_eV = 23.33;
    static constexpr double SIII_eV = 34.83;
    static constexpr double SIV_eV = 47.31;

    static constexpr double FeI_eV = 7.902;
    static constexpr double FeII_eV = 16.19;
    static constexpr double FeIII_eV = 30.65;
    static constexpr double FeIV_eV = 54.8;
    static constexpr double FeV_eV = 75.01;
    static constexpr double FeVI_eV = 99.06;
    static constexpr double FeVII_eV = 125.;
    static constexpr double FeVIII_eV = 151.1;
    static constexpr double FeIX_eV = 233.6;
    static constexpr double FeX_eV = 262.1;
    static constexpr double FeXI_eV = 290.2;
    static constexpr double FeXII_eV = 330.8;
    static constexpr double FeXIII_eV = 361.;
    static constexpr double FeXIV_eV = 392.2;
    static constexpr double FeXV_eV = 457.;
    static constexpr double FeXVI_eV = 489.3;

    // total number of ions with cross-sections
    static constexpr int numIons = 74;

    // ---- element/ion index ranges ----

    static constexpr int firstH = 0;
    static constexpr int numH = 1;  // HI only (HII is fully ionized, no XS)
    static constexpr int firstHe = 1;
    static constexpr int numHe = 2;  // HeI, HeII
    static constexpr int firstC = 3;
    static constexpr int numC = 6;  // CI-CVI
    static constexpr int firstN = 9;
    static constexpr int numN = 7;  // NI-NVII
    static constexpr int firstO = 16;
    static constexpr int numO = 8;  // OI-OVIII
    static constexpr int firstNe = 24;
    static constexpr int numNe = 8;  // NeI-NeVIII
    static constexpr int firstMg = 32;
    static constexpr int numMg = 10;  // MgI-MgX
    static constexpr int firstSi = 42;
    static constexpr int numSi = 12;  // SiI-SiXII
    static constexpr int firstS = 54;
    static constexpr int numS = 4;  // SI-SIV
    static constexpr int firstFe = 58;
    static constexpr int numFe = 16;  // FeI-FeXVI

    // ---- dispatch by linear ion index ----

    /** Ion ordering convention used throughout the photoionization solver.
        Index 0-2:   H I, He I, He II
        Index 3-8:   C I-VI
        Index 9-15:  N I-VII
        Index 16-23: O I-VIII
        Index 24-31: Ne I-VIII
        Index 32-41: Mg I-X
        Index 42-53: Si I-XII
        Index 54-57: S I-IV
        Index 58-73: Fe I-XVI */

    /** Returns the ionization potential [eV] for the given linear ion index. */
    static double ionizationPotential(int ion);

    /** Returns the photoionization cross-section [cm^2] for the given linear ion index
        at photon energy E [eV]. */
    static double sigma(int ion, double E);

    /** HI photoionization cross-section [cm^2] over [13.6, 50000] eV. */
    static double sigmaHI(double E);

    /** HeI photoionization cross-section [cm^2] over [24.59, 50000] eV. */
    static double sigmaHeI(double E);

    /** HeII photoionization cross-section [cm^2] over [54.42, 50000] eV. */
    static double sigmaHeII(double E);
};

//////////////////////////////////////////////////////////////////////

#endif
