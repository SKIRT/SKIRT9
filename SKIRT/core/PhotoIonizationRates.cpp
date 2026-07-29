/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PhotoIonizationRates.hpp"
#include "VernerCrossSections.hpp"

//////////////////////////////////////////////////////////////////////

namespace
{
    // conversion factor: eV -> K (IP_eV * eV_to_K = IP in Kelvin)
    constexpr double eV_to_K = 11604.52;

    // ============================================================
    //  Hydrogen
    // ============================================================

    /** Case B recombination coefficient for HII [cm^3/s]. Hui & Gnedin (1997). */
    inline double alphaBHII(double T)
    {
        double lambda = 315614. / T;
        return 2.753e-14 * pow(lambda, 1.5) / pow(1. + pow(lambda / 2.74, 0.407), 2.242);
    }

    /** HI collisional ionization rate [cm^3/s]. Voronov (1997). */
    inline double betaHI(double T)
    {
        double U = 13.598 * eV_to_K / T;
        return 2.91e-8 / (0.232 + U) * pow(U, 0.39) * exp(-U);
    }

    // ============================================================
    //  Helium
    // ============================================================

    /** Case B recombination coefficient for HeII [cm^3/s]. */
    inline double alphaBHeII(double T)
    {
        double lambda = 570670. / T;
        return 1.26e-14 * pow(lambda, 0.75);
    }

    /** Dielectronic recombination coefficient for HeII [cm^3/s]. */
    inline double alphaDRHeII(double T)
    {
        return 1.9e-3 * pow(T, -1.5) * exp(-470000. / T) * (1. + 0.3 * exp(-94000. / T));
    }

    /** Case B recombination coefficient for HeIII [cm^3/s]. Hui & Gnedin (1997). */
    inline double alphaBHeIII(double T)
    {
        double lambda = 1263030. / T;
        return 5.506e-14 * pow(lambda, 1.5) / pow(1. + pow(lambda / 2.74, 0.407), 2.242);
    }

    /** HeI collisional ionization rate [cm^3/s]. Voronov (1997). */
    inline double betaHeI(double T)
    {
        double U = 24.587 * eV_to_K / T;
        return 9.38e-9 / (0.231 + U) * pow(U, 0.39) * exp(-U);
    }

    /** HeII collisional ionization rate [cm^3/s]. Voronov (1997). */
    inline double betaHeII(double T)
    {
        double U = 54.418 * eV_to_K / T;
        return 5.68e-8 / (0.265 + U) * pow(U, 0.25) * exp(-U);
    }

    /** HeIII + HI charge exchange recombination [cm^3/s]. */
    inline double xiHIHeIII(double /*T*/)
    {
        return 1e-14;
    }

    /** HeII + HI charge exchange recombination [cm^3/s]. */
    inline double xiHIHeII(double T)
    {
        double T4 = (T < 6e3) ? 0.6 : ((T > 1e5) ? 10. : 1e-4 * T);
        return 7.47e-15 * pow(T4, 2.06) * (1. + 9.93 * exp(-3.89 * T4));
    }

    // ============================================================
    //  Carbon  (CI-CVI: 6 ionization stages -> CVII is bare nucleus)
    // ============================================================

    // -- radiative recombination (ion+1 -> ion) --

    inline double alphaRRCVII(double T)
    {
        double s0 = sqrt(T / 95.02), s1 = sqrt(T / 2.517e7);
        return 5.337e-10 / (s0 * pow(1. + s0, 0.2515) * pow(1. + s1, 1.7485));
    }
    inline double alphaRRCVI(double T)
    {
        double s0 = sqrt(T / 264.7), s1 = sqrt(T / 2.773e7);
        return 2.044e-10 / (s0 * pow(1. + s0, 0.3258) * pow(1. + s1, 1.6742));
    }
    inline double alphaRRCV(double T)
    {
        double s0 = sqrt(T / 1355.), s1 = sqrt(T / 1.872e7);
        return 4.798e-11 / (s0 * pow(1. + s0, 0.5166) * pow(1. + s1, 1.4834));
    }
    inline double alphaRRCIV(double T)
    {
        double s0 = sqrt(T / 111.5), s1 = sqrt(T / 5.938e6);
        return 1.12e-10 / (s0 * pow(1. + s0, 0.3263) * pow(1. + s1, 1.6737));
    }
    inline double alphaRRCIII(double T)
    {
        double s0 = sqrt(T / 0.1643), s1 = sqrt(T / 2.172e6);
        double Bp = 0.8012 + 0.0427 * exp(-63410. / T);
        return 2.067e-9 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRCII(double T)
    {
        double s0 = sqrt(T / 0.00667), s1 = sqrt(T / 1.943e6);
        double Bp = 0.7849 + 0.1597 * exp(-49550. / T);
        return 2.995e-9 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }

    // -- dielectronic recombination --

    inline double alphaDRCVI(double T)
    {
        return pow(T, -1.5)
               * (0.001426 * exp(-3.116e6 / T) + 0.03046 * exp(-4.075e6 / T) + 0.0008373 * exp(-5.749e6 / T));
    }
    inline double alphaDRCV(double T)
    {
        return pow(T, -1.5)
               * (0.002646 * exp(-2.804e6 / T) + 0.01762 * exp(-3.485e6 / T) - 0.0007843 * exp(-4.324e6 / T));
    }
    inline double alphaDRCIV(double T)
    {
        return pow(T, -1.5)
               * (4.673e-7 * exp(-723.3 / T) + 1.887e-5 * exp(-2847. / T) + 1.305e-5 * exp(-10540. / T)
                  + 0.003099 * exp(-89150. / T) + 0.0003001 * exp(-281200. / T) + 0.002553 * exp(-3.254e6 / T));
    }
    inline double alphaDRCIII(double T)
    {
        return pow(T, -1.5)
               * (3.489e-6 * exp(-2660. / T) + 2.222e-7 * exp(-3756. / T) + 1.954e-5 * exp(-25660. / T)
                  + 0.004212 * exp(-140000. / T) + 0.0002037 * exp(-1.801e6 / T) + 0.0002936 * exp(-4.307e6 / T));
    }
    inline double alphaDRCII(double T)
    {
        return pow(T, -1.5)
               * (6.346e-9 * exp(-12.17 / T) + 9.793e-9 * exp(-73.8 / T) + 1.634e-6 * exp(-15230. / T)
                  + 0.0008369 * exp(-120700. / T) + 0.0003355 * exp(-214400. / T));
    }

    // -- collisional ionization (Voronov 1997) --

    inline double betaCI(double T)
    {
        double U = VernerCrossSections::CI_eV * eV_to_K / T;
        return 6.85e-8 * pow(U, 0.25) * exp(-U) / (0.193 + U);
    }
    inline double betaCII(double T)
    {
        double U = VernerCrossSections::CII_eV * eV_to_K / T;
        return 1.86e-8 * (1. + sqrt(U)) * pow(U, 0.24) * exp(-U) / (0.286 + U);
    }
    inline double betaCIII(double T)
    {
        double U = VernerCrossSections::CIII_eV * eV_to_K / T;
        return 6.35e-9 * (1. + sqrt(U)) * pow(U, 0.21) * exp(-U) / (0.427 + U);
    }
    inline double betaCIV(double T)
    {
        double U = VernerCrossSections::CIV_eV * eV_to_K / T;
        return 1.5e-9 * (1. + sqrt(U)) * pow(U, 0.13) * exp(-U) / (0.416 + U);
    }
    inline double betaCV(double T)
    {
        double U = VernerCrossSections::CV_eV * eV_to_K / T;
        return 2.99e-10 * (1. + sqrt(U)) * pow(U, 0.02) * exp(-U) / (0.666 + U);
    }
    inline double betaCVI(double T)
    {
        double U = VernerCrossSections::CVI_eV * eV_to_K / T;
        return 1.23e-10 * (1. + sqrt(U)) * pow(U, 0.16) * exp(-U) / (0.62 + U);
    }

    // -- charge exchange --

    inline double xiHICV(double T)
    {
        double T4 = (T < 10.) ? 0.001 : ((T > 1e5) ? 10. : 1e-4 * T);
        return 3.3246e-7 * pow(T4, -0.11) * (1. - 0.995 * exp(-0.00158 * T4));
    }
    inline double xiHICIV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 1e5) ? 10. : 1e-4 * T);
        return 3.25e-9 * pow(T4, 0.21) * (1. + 0.19 * exp(-3.29 * T4));
    }
    inline double xiHICIII(double T)
    {
        double T4 = (T < 5e3) ? 0.5 : ((T > 5e4) ? 5. : 1e-4 * T);
        return 1.67e-13 * pow(T4, 2.79) * (1. + 304.72 * exp(-4.07 * T4));
    }
    inline double xiHICII(double T)
    {
        double T4 = (T < 5500.) ? 0.55 : ((T > 1e5) ? 10. : 1e-4 * T);
        return 4.88e-16 * pow(T4, 3.25) * (1. - 1.12 * exp(-0.21 * T4));
    }
    inline double chiHIICI(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 1e5) ? 10. : 1e-4 * T);
        return 1.07e-15 * pow(T4, 3.15) * (1. + 176.43 * exp(-4.29 * T4));
    }

    // ============================================================
    //  Nitrogen  (NI-NVII: 7 ionization stages -> NVIII is bare nucleus)
    // ============================================================

    inline double alphaRRNVIII(double T)
    {
        double s0 = sqrt(T / 131.6), s1 = sqrt(T / 3.427e7);
        return 6.17e-10 / (s0 * pow(1. + s0, 0.2519) * pow(1. + s1, 1.7481));
    }
    inline double alphaRRNVII(double T)
    {
        double s0 = sqrt(T / 396.), s1 = sqrt(T / 3.583e7);
        return 2.388e-10 / (s0 * pow(1. + s0, 0.3268) * pow(1. + s1, 1.6732));
    }
    inline double alphaRRNVI(double T)
    {
        double s0 = sqrt(T / 1957.), s1 = sqrt(T / 2.177e7);
        return 6.245e-11 / (s0 * pow(1. + s0, 0.5015) * pow(1. + s1, 1.4985));
    }
    inline double alphaRRNV(double T)
    {
        double s0 = sqrt(T / 182.3), s1 = sqrt(T / 7.751e6);
        return 1.533e-10 / (s0 * pow(1. + s0, 0.3318) * pow(1. + s1, 1.6682));
    }
    inline double alphaRRNIV(double T)
    {
        double s0 = sqrt(T / 3.75), s1 = sqrt(T / 3.468e6);
        double Bp = 0.7768 + 0.0223 * exp(-72060. / T);
        return 7.923e-10 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRNIII(double T)
    {
        double s0 = sqrt(T / 0.1231), s1 = sqrt(T / 3.016e6);
        double Bp = 0.7948 + 0.0774 * exp(-101600. / T);
        return 2.41e-9 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRNII(double T)
    {
        double s0 = sqrt(T / 0.09467), s1 = sqrt(T / 2.954e6);
        double Bp = 0.7308 + 0.244 * exp(-67390. / T);
        return 6.387e-10 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }

    inline double alphaDRNVII(double T)
    {
        return pow(T, -1.5)
               * (0.002801 * exp(-4.198e6 / T) + 0.04362 * exp(-5.516e6 / T) + 0.001117 * exp(-8.05e6 / T));
    }
    inline double alphaDRNVI(double T)
    {
        return pow(T, -1.5) * (0.005761 * exp(-3.86e6 / T) + 0.03434 * exp(-4.883e6 / T) - 0.00166 * exp(-6.259e6 / T));
    }
    inline double alphaDRNV(double T)
    {
        return pow(T, -1.5)
               * (2.04e-6 * exp(-3084. / T) + 6.986e-5 * exp(-13320. / T) + 0.0003168 * exp(-64750. / T)
                  + 0.004353 * exp(-118100. / T) + 0.0007765 * exp(-668700. / T) + 0.005101 * exp(-4.778e6 / T));
    }
    inline double alphaDRNIV(double T)
    {
        return pow(T, -1.5)
               * (3.386e-6 * exp(-1406. / T) + 3.036e-5 * exp(-6965. / T) + 5.945e-5 * exp(-26040. / T)
                  + 0.001195 * exp(-130400. / T) + 0.006462 * exp(-196500. / T) + 0.001358 * exp(-4.466e6 / T));
    }
    inline double alphaDRNIII(double T)
    {
        return pow(T, -1.5)
               * (7.712e-8 * exp(-71.13 / T) + 4.839e-8 * exp(-276.5 / T) + 2.218e-6 * exp(-14390. / T)
                  + 0.001536 * exp(-134700. / T) + 0.003647 * exp(-249600. / T) + 4.234e-5 * exp(-2.204e6 / T));
    }
    inline double alphaDRNII(double T)
    {
        return pow(T, -1.5)
               * (1.658e-8 * exp(-12.65 / T) + 2.76e-8 * exp(-84.25 / T) + 2.391e-9 * exp(-296.4 / T)
                  + 7.585e-7 * exp(-5923. / T) + 0.0003012 * exp(-127800. / T) + 0.0007132 * exp(-218400. / T));
    }

    inline double betaNI(double T)
    {
        double U = VernerCrossSections::NI_eV * eV_to_K / T;
        return 4.82e-8 * pow(U, 0.42) * exp(-U) / (0.0652 + U);
    }
    inline double betaNII(double T)
    {
        double U = VernerCrossSections::NII_eV * eV_to_K / T;
        return 2.98e-8 * pow(U, 0.3) * exp(-U) / (0.31 + U);
    }
    inline double betaNIII(double T)
    {
        double U = VernerCrossSections::NIII_eV * eV_to_K / T;
        return 8.1e-9 * (1. + sqrt(U)) * pow(U, 0.24) * exp(-U) / (0.35 + U);
    }
    inline double betaNIV(double T)
    {
        double U = VernerCrossSections::NIV_eV * eV_to_K / T;
        return 3.71e-9 * (1. + sqrt(U)) * pow(U, 0.18) * exp(-U) / (0.549 + U);
    }
    inline double betaNV(double T)
    {
        double U = VernerCrossSections::NV_eV * eV_to_K / T;
        return 1.51e-9 * pow(U, 0.74) * exp(-U) / (0.0167 + U);
    }
    inline double betaNVI(double T)
    {
        double U = VernerCrossSections::NVI_eV * eV_to_K / T;
        return 3.71e-10 * pow(U, 0.29) * exp(-U) / (0.546 + U);
    }
    inline double betaNVII(double T)
    {
        double U = VernerCrossSections::NVII_eV * eV_to_K / T;
        return 7.77e-11 * (1. + sqrt(U)) * pow(U, 0.16) * exp(-U) / (0.624 + U);
    }

    inline double xiHINV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 1e6) ? 100. : 1e-4 * T);
        return 2.95e-9 * pow(T4, 0.55) * (1. - 0.39 * exp(-1.07 * T4));
    }
    inline double xiHINIV(double T)
    {
        double T4 = (T < 10.) ? 0.001 : ((T > 1e5) ? 10. : 1e-4 * T);
        return 4.54e-9 * pow(T4, 0.57) * (1. - 0.65 * exp(-0.89 * T4));
    }
    inline double xiHINIII(double T)
    {
        if (T < 1500.) return 0.8692e-9 * pow(T / 1500., 0.17);
        if (T < 20000.) return 0.9703e-9 * pow(T / 10000., 0.058);
        return 1.0101e-9 + 1.4589e-9 * pow(log10(T / 20000.), 2.06);
    }
    inline double xiHINII(double T)
    {
        double lnT = log(T);
        return ((1.56e-15 - 1.79e-16 * lnT + 1.15e-20 * T) * T + 1.08e-13 * lnT)
               + ((6.83e-16 - 7.40e-17 * lnT + 3.73e-21 * T) * T + 1.75e-15 * lnT) * exp(-16680. / T);
    }
    inline double chiHIINI(double T)
    {
        double lnT = log(T);
        return ((1.64e-16 - 8.76e-17 * lnT + 2.41e-20 * T) * T + 9.83e-13 * lnT) * exp(-10985. / T);
    }

    // ============================================================
    //  Oxygen  (OI-OVIII: 8 ionization stages -> OIX is bare nucleus)
    // ============================================================

    inline double alphaRROIX(double T)
    {
        double s0 = sqrt(T / 195.1), s1 = sqrt(T / 4.483e7);
        return 6.552e-10 / (s0 * pow(1. + s0, 0.253) * pow(1. + s1, 1.747));
    }
    inline double alphaRROVIII(double T)
    {
        double s0 = sqrt(T / 584.2), s1 = sqrt(T / 4.559e7);
        return 2.652e-10 / (s0 * pow(1. + s0, 0.3295) * pow(1. + s1, 1.6705));
    }
    inline double alphaRROVII(double T)
    {
        double s0 = sqrt(T / 2392.), s1 = sqrt(T / 2.487e7);
        return 8.193e-11 / (s0 * pow(1. + s0, 0.4835) * pow(1. + s1, 1.5165));
    }
    inline double alphaRROVI(double T)
    {
        double s0 = sqrt(T / 337.2), s1 = sqrt(T / 1.03e7);
        return 1.724e-10 / (s0 * pow(1. + s0, 0.3444) * pow(1. + s1, 1.6556));
    }
    inline double alphaRROV(double T)
    {
        double s0 = sqrt(T / 0.6821), s1 = sqrt(T / 5.076e6);
        return 3.955e-9 / (s0 * pow(1. + s0, 0.2187) * pow(1. + s1, 1.7813));
    }
    inline double alphaRROIV(double T)
    {
        double s0 = sqrt(T / 0.5235), s1 = sqrt(T / 4.47e6);
        double Bp = 0.7844 + 0.0447 * exp(-164200. / T);
        return 2.501e-9 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRROIII(double T)
    {
        double s0 = sqrt(T / 0.1602), s1 = sqrt(T / 4.377e6);
        double Bp = 0.7668 + 0.107 * exp(-139200. / T);
        return 2.096e-9 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRROII(double T)
    {
        double s0 = sqrt(T / 4.136), s1 = sqrt(T / 4.214e6);
        double Bp = 0.6109 + 0.4093 * exp(-87700. / T);
        return 6.622e-11 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }

    inline double alphaDROVIII(double T)
    {
        return pow(T, -1.5) * (0.004925 * exp(-5.44e6 / T) + 0.05837 * exp(-7.17e6 / T) + 0.001359 * exp(-1.152e7 / T));
    }
    inline double alphaDROVII(double T)
    {
        return pow(T, -1.5) * (0.06135 * exp(-6.113e6 / T) + 0.0001968 * exp(-3.656e7 / T));
    }
    inline double alphaDROVI(double T)
    {
        return pow(T, -1.5)
               * (2.389e-5 * exp(-23260. / T) + 0.0001355 * exp(-32090. / T) + 0.005885 * exp(-131600. / T)
                  + 0.002163 * exp(-673100. / T) + 0.0006341 * exp(-1.892e6 / T) + 0.01348 * exp(-6.15e6 / T));
    }
    inline double alphaDROV(double T)
    {
        return pow(T, -1.5)
               * (1.615e-5 * exp(-756.9 / T) + 9.299e-6 * exp(-3659. / T) + 0.000153 * exp(-19840. / T)
                  + 0.0006616 * exp(-84290. / T) + 0.0108 * exp(-229400. / T) + 0.0007503 * exp(-1.161e6 / T)
                  + 0.002892 * exp(-6.137e6 / T));
    }
    inline double alphaDROIV(double T)
    {
        return pow(T, -1.5)
               * (3.932e-7 * exp(-150.9 / T) + 2.523e-7 * exp(-621.1 / T) + 3.447e-5 * exp(-15620. / T)
                  + 0.005776 * exp(-193600. / T) + 0.005101 * exp(-470000. / T));
    }
    inline double alphaDROIII(double T)
    {
        return pow(T, -1.5)
               * (1.627e-7 * exp(-45.35 / T) + 1.262e-7 * exp(-284.7 / T) + 6.663e-7 * exp(-4166. / T)
                  + 3.925e-6 * exp(-28770. / T) + 0.002406 * exp(-195300. / T) + 0.001146 * exp(-364600. / T));
    }
    inline double alphaDROII(double T)
    {
        return pow(T, -1.5)
               * (5.629e-8 * exp(-5395. / T) + 2.55e-7 * exp(-17700. / T) + 0.0006173 * exp(-167100. / T)
                  + 0.0001627 * exp(-268700. / T));
    }

    inline double betaOI(double T)
    {
        double U = VernerCrossSections::OI_eV * eV_to_K / T;
        return 3.59e-8 * pow(U, 0.34) * exp(-U) / (0.073 + U);
    }
    inline double betaOII(double T)
    {
        double U = VernerCrossSections::OII_eV * eV_to_K / T;
        return 1.39e-8 * (1. + sqrt(U)) * pow(U, 0.22) * exp(-U) / (0.212 + U);
    }
    inline double betaOIII(double T)
    {
        double U = VernerCrossSections::OIII_eV * eV_to_K / T;
        return 9.31e-9 * (1. + sqrt(U)) * pow(U, 0.27) * exp(-U) / (0.27 + U);
    }
    inline double betaOIV(double T)
    {
        double U = VernerCrossSections::OIV_eV * eV_to_K / T;
        return 1.02e-8 * pow(U, 0.27) * exp(-U) / (0.614 + U);
    }
    inline double betaOV(double T)
    {
        double U = VernerCrossSections::OV_eV * eV_to_K / T;
        return 2.19e-9 * (1. + sqrt(U)) * pow(U, 0.17) * exp(-U) / (0.63 + U);
    }
    inline double betaOVI(double T)
    {
        double U = VernerCrossSections::OVI_eV * eV_to_K / T;
        return 1.95e-9 * pow(U, 0.54) * exp(-U) / (0.36 + U);
    }
    inline double betaOVII(double T)
    {
        double U = VernerCrossSections::OVII_eV * eV_to_K / T;
        return 2.12e-10 * pow(U, 0.35) * exp(-U) / (0.396 + U);
    }
    inline double betaOVIII(double T)
    {
        double U = VernerCrossSections::OVIII_eV * eV_to_K / T;
        return 5.21e-11 * (1. + sqrt(U)) * pow(U, 0.16) * exp(-U) / (0.629 + U);
    }

    inline double xiHIOV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 3e4) ? 3. : 1e-4 * T);
        return 2.52e-10 * pow(T4, 0.63) * (1. + 2.08 * exp(-4.16 * T4));
    }
    inline double xiHIOIV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 5e4) ? 5. : 1e-4 * T);
        return 3.98e-9 * pow(T4, 0.26) * (1. + 0.56 * exp(-2.62 * T4));
    }
    inline double xiHIOIII(double T)
    {
        if (T < 1500.) return 0.5337e-9 * pow(T / 100., -0.076);
        return 0.4344e-9 + 0.6340e-9 * pow(log10(T / 1500.), 2.06);
    }
    inline double xiHIOII(double T)
    {
        if (T < 10.) return 3.744e-10;
        double lnT = log(T);
        return ((((1.1963502e-13 * lnT - 2.8577012e-12) * lnT + 2.9979994e-11) * lnT - 1.3146803e-10) * lnT
                + 2.3651505e-10)
                   * lnT
               + 2.3344302e-10;
    }
    inline double chiHIIOI(double T)
    {
        if (T < 10.) return 4.749e-20;
        if (T < 190.) return exp(-21.134531 - 242.06831 / T + 84.761441 / (T * T));
        if (T < 200.) return 2.18733e-12 * (T - 190.) + 1.85823e-10;
        double lnT = log(T);
        return (((((1.1580844e-14 * lnT - 2.6139493e-13) * lnT + 2.0699463e-12) * lnT - 3.6606214e-12) * lnT
                 - 1.488594e-12)
                    * lnT
                - 3.7282001e-13)
                   * lnT
               - 7.6767404e-14;
    }

    // ============================================================
    //  Neon  (NeI-NeVIII: 8 ionization stages -> NeIX is bare nucleus)
    // ============================================================

    inline double alphaRRNeIX(double T)
    {
        double s0 = sqrt(T / 3647.), s1 = sqrt(T / 3.365e7);
        return 1.186e-10 / (s0 * pow(1. + s0, 0.4646) * pow(1. + s1, 1.5354));
    }
    inline double alphaRRNeVIII(double T)
    {
        double s0 = sqrt(T / 510.2), s1 = sqrt(T / 1.535e7);
        return 2.755e-10 / (s0 * pow(1. + s0, 0.3414) * pow(1. + s1, 1.6586));
    }
    inline double alphaRRNeVII(double T)
    {
        double s0 = sqrt(T / 6.293), s1 = sqrt(T / 8.091e6);
        return 2.557e-9 / (s0 * pow(1. + s0, 0.2399) * pow(1. + s1, 1.7601));
    }
    inline double alphaRRNeVI(double T)
    {
        double s0 = sqrt(T / 13.11), s1 = sqrt(T / 8.047e6);
        double Bp = 0.7556 + 0.025 * exp(-277100. / T);
        return 1.127e-9 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRNeV(double T)
    {
        double s0 = sqrt(T / 2.504), s1 = sqrt(T / 8.037e6);
        double Bp = 0.7593 + 0.0406 * exp(-325500. / T);
        return 1.861e-9 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRNeIV(double T)
    {
        double s0 = sqrt(T / 3.332), s1 = sqrt(T / 8.696e6);
        double Bp = 0.7254 + 0.0921 * exp(-304400. / T);
        return 8.321e-10 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRNeIII(double T)
    {
        double s0 = sqrt(T / 9.924), s1 = sqrt(T / 8.878e6);
        double Bp = 0.6434 + 0.2205 * exp(-229200. / T);
        return 1.773e-10 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRNeII(double T)
    {
        double s0 = sqrt(T / 67.39), s1 = sqrt(T / 7.563e6);
        double Bp = 0.3556 + 0.6472 * exp(-159800. / T);
        return 1.295e-11 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }

    inline double alphaDRNeVIII(double T)
    {
        return pow(T, -1.5)
               * (0.0002399 * exp(-25360. / T) + 0.0003532 * exp(-48000. / T) + 0.008928 * exp(-177000. / T)
                  + 0.005427 * exp(-1.03e6 / T) + 0.005342 * exp(-1.859e6 / T) + 0.03981 * exp(-9.743e6 / T));
    }
    inline double alphaDRNeVII(double T)
    {
        return pow(T, -1.5)
               * (8.46e-5 * exp(-1049. / T) + 0.0001817 * exp(-3829. / T) + 0.001176 * exp(-61330. / T)
                  + 0.01397 * exp(-256800. / T) + 0.005566 * exp(-460000. / T) + 0.006229 * exp(-1.324e6 / T)
                  + 0.01031 * exp(-9.353e6 / T));
    }
    inline double alphaDRNeVI(double T)
    {
        return pow(T, -1.5)
               * (5.653e-6 * exp(-628. / T) + 4.344e-5 * exp(-2812. / T) + 0.0001086 * exp(-13240. / T)
                  + 0.000598 * exp(-80640. / T) + 0.01457 * exp(-305200. / T) + 0.01601 * exp(-1.032e6 / T)
                  + 0.0005365 * exp(-2.388e6 / T));
    }
    inline double alphaDRNeV(double T)
    {
        return pow(T, -1.5)
               * (2.922e-6 * exp(-205. / T) + 7.144e-6 * exp(-2205. / T) + 2.836e-5 * exp(-9271. / T)
                  + 9.82e-5 * exp(-49880. / T) + 0.008379 * exp(-290400. / T) + 0.01009 * exp(-878200. / T));
    }
    inline double alphaDRNeIV(double T)
    {
        return pow(T, -1.5)
               * (2.763e-6 * exp(-639.3 / T) + 1.053e-5 * exp(-1499. / T) + 4.453e-5 * exp(-32270. / T)
                  + 0.006244 * exp(-256100. / T) + 0.0003146 * exp(-450500. / T) + 0.004465 * exp(-793400. / T));
    }
    inline double alphaDRNeIII(double T)
    {
        return pow(T, -1.5)
               * (2.98e-8 * exp(-45.79 / T) + 1.257e-7 * exp(-475.3 / T) + 1.122e-6 * exp(-14810. / T)
                  + 0.002626 * exp(-281000. / T) + 0.0008802 * exp(-476300. / T) + 1.231 * exp(-4.677e8 / T));
    }
    inline double alphaDRNeII(double T)
    {
        return pow(T, -1.5)
               * (4.152e-9 * exp(-26.89 / T) + 4.656e-9 * exp(-202.1 / T) + 1.31e-8 * exp(-720. / T)
                  + 1.417e-9 * exp(-48920. / T) + 0.0007968 * exp(-314400. / T) + 1.271e-5 * exp(-673800. / T));
    }

    inline double betaNeI(double T)
    {
        double U = VernerCrossSections::NeI_eV * eV_to_K / T;
        return 1.5e-8 * (1. + sqrt(U)) * pow(U, 0.43) * exp(-U) / (0.0329 + U);
    }
    inline double betaNeII(double T)
    {
        double U = VernerCrossSections::NeII_eV * eV_to_K / T;
        return 1.98e-8 * pow(U, 0.2) * exp(-U) / (0.295 + U);
    }
    inline double betaNeIII(double T)
    {
        double U = VernerCrossSections::NeIII_eV * eV_to_K / T;
        return 7.03e-9 * (1. + sqrt(U)) * pow(U, 0.39) * exp(-U) / (0.0677 + U);
    }
    inline double betaNeIV(double T)
    {
        double U = VernerCrossSections::NeIV_eV * eV_to_K / T;
        return 4.24e-9 * (1. + sqrt(U)) * pow(U, 0.58) * exp(-U) / (0.0482 + U);
    }
    inline double betaNeV(double T)
    {
        double U = VernerCrossSections::NeV_eV * eV_to_K / T;
        return 2.79e-9 * (1. + sqrt(U)) * pow(U, 0.25) * exp(-U) / (0.305 + U);
    }
    inline double betaNeVI(double T)
    {
        double U = VernerCrossSections::NeVI_eV * eV_to_K / T;
        return 3.45e-9 * pow(U, 0.28) * exp(-U) / (0.581 + U);
    }
    inline double betaNeVII(double T)
    {
        double U = VernerCrossSections::NeVII_eV * eV_to_K / T;
        return 9.56e-10 * (1. + sqrt(U)) * pow(U, 0.14) * exp(-U) / (0.749 + U);
    }
    inline double betaNeVIII(double T)
    {
        double U = VernerCrossSections::NeVIII_eV * eV_to_K / T;
        return 4.73e-10 * (1. + sqrt(U)) * pow(U, 0.04) * exp(-U) / (0.992 + U);
    }

    inline double xiHINeV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 3e4) ? 3. : 1e-4 * T);
        return 6.47e-9 * pow(T4, 0.54) * (1. + 3.59 * exp(-5.22 * T4));
    }
    inline double xiHINeIV(double T)
    {
        double T4 = (T < 5e3) ? 0.5 : ((T > 5e4) ? 5. : 1e-4 * T);
        return 1.473e-8 * pow(T4, 0.0452) * (1. - 0.84 * exp(-0.31 * T4));
    }
    inline double xiHINeIII(double /*T*/)
    {
        return 1e-14;
    }

    // ============================================================
    //  Magnesium  (MgI-MgX: 10 ionization stages, rates up to MgXI)
    // ============================================================

    // -- radiative recombination (ion+1 -> ion) --

    inline double alphaRRMgXI(double T)
    {
        double s0 = sqrt(T / 4944.), s1 = sqrt(T / 4.434e7);
        return 1.602e-10 / (s0 * pow(1. + s0, 0.4508) * pow(1. + s1, 1.5492));
    }
    inline double alphaRRMgX(double T)
    {
        double s0 = sqrt(T / 869.3), s1 = sqrt(T / 2.196e7);
        return 3.445e-10 / (s0 * pow(1. + s0, 0.3447) * pow(1. + s1, 1.6553));
    }
    inline double alphaRRMgIX(double T)
    {
        double s0 = sqrt(T / 38.06), s1 = sqrt(T / 1.205e7);
        return 1.634e-9 / (s0 * pow(1. + s0, 0.2565) * pow(1. + s1, 1.7435));
    }
    inline double alphaRRMgVIII(double T)
    {
        double s0 = sqrt(T / 5.587), s1 = sqrt(T / 1.235e7);
        return 3.859e-9 / (s0 * pow(1. + s0, 0.2421) * pow(1. + s1, 1.7579));
    }
    inline double alphaRRMgVII(double T)
    {
        double s0 = sqrt(T / 36.74), s1 = sqrt(T / 1.263e7);
        double Bp = 0.7353 + 0.0211 * exp(-404900. / T);
        return 9.133e-10 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRMgVI(double T)
    {
        double s0 = sqrt(T / 25.82), s1 = sqrt(T / 1.355e7);
        double Bp = 0.7203 + 0.0436 * exp(-569100. / T);
        return 7.515e-10 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRMgV(double T)
    {
        double s0 = sqrt(T / 32.05), s1 = sqrt(T / 1.626e7);
        double Bp = 0.6803 + 0.0764 * exp(-539900. / T);
        return 4.031e-10 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRMgIV(double T)
    {
        double s0 = sqrt(T / 77.48), s1 = sqrt(T / 2.015e7);
        double Bp = 0.56 + 0.1917 * exp(-513900. / T);
        return 1.249e-10 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRMgIII(double T)
    {
        double s0 = sqrt(T / 787.7), s1 = sqrt(T / 7.925e7);
        double Bp = 0.1074 + 0.4631 * exp(-502700. / T);
        return 1.345e-11 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRMgII(double T)
    {
        double s0 = sqrt(T / 5.637), s1 = sqrt(T / 1.551e6);
        double Bp = 0.6845 + 0.3945 * exp(-836000. / T);
        return 5.452e-11 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }

    // -- dielectronic recombination --

    inline double alphaDRMgXI(double T)
    {
        return pow(T, -1.5) * (0.03067 * exp(-1.136e7 / T) + 0.1375 * exp(-1.431e7 / T) + 0.01347 * exp(-1.762e7 / T));
    }
    inline double alphaDRMgX(double T)
    {
        return pow(T, -1.5)
               * (0.0002582 * exp(-28610. / T) + 0.0005033 * exp(-40870. / T) + 0.004561 * exp(-152800. / T)
                  + 0.008754 * exp(-267700. / T) + 0.02734 * exp(-2.055e6 / T) + 0.07509 * exp(-1.298e7 / T)
                  + 0.01385 * exp(-1.852e7 / T));
    }
    inline double alphaDRMgIX(double T)
    {
        return pow(T, -1.5)
               * (3.565e-5 * exp(-1759. / T) + 2.017e-5 * exp(-12370. / T) + 0.001165 * exp(-56510. / T)
                  + 0.005016 * exp(-168000. / T) + 0.02434 * exp(-400700. / T) + 0.02508 * exp(-1.991e6 / T)
                  + 0.02475 * exp(-1.33e7 / T) + 0.0007087 * exp(-4.06e7 / T));
    }
    inline double alphaDRMgVIII(double T)
    {
        return pow(T, -1.5)
               * (5.451e-5 * exp(-219.6 / T) + 6.999e-5 * exp(-4524. / T) + 0.0003928 * exp(-24910. / T)
                  + 0.0007483 * exp(-91110. / T) + 0.01851 * exp(-336700. / T) + 0.0119 * exp(-831600. / T)
                  + 0.04764 * exp(-1.817e6 / T));
    }
    inline double alphaDRMgVII(double T)
    {
        return pow(T, -1.5)
               * (2.931e-5 * exp(-766.7 / T) + 5.192e-5 * exp(-3444. / T) + 0.0001769 * exp(-24090. / T)
                  + 0.0007988 * exp(-119000. / T) + 0.01488 * exp(-399900. / T) + 0.04546 * exp(-1.502e6 / T)
                  + 0.0001505 * exp(-2.464e7 / T));
    }
    inline double alphaDRMgVI(double T)
    {
        return pow(T, -1.5)
               * (2.407e-5 * exp(-8061. / T) + 0.0001187 * exp(-21540. / T) + 0.0003845 * exp(-66400. / T)
                  + 0.01333 * exp(-340600. / T) + 0.02804 * exp(-1.327e6 / T) + 0.0008746 * exp(-3.347e6 / T));
    }
    inline double alphaDRMgV(double T)
    {
        return pow(T, -1.5)
               * (2.574e-6 * exp(-112.9 / T) + 2.001e-6 * exp(-1585. / T) + 1.829e-5 * exp(-11930. / T)
                  + 2.362e-5 * exp(-41510. / T) + 0.006413 * exp(-368400. / T) + 0.004118 * exp(-692500. / T)
                  + 0.007224 * exp(-1.15e6 / T));
    }
    inline double alphaDRMgIV(double T)
    {
        return pow(T, -1.5)
               * (9.4e-8 * exp(-266. / T) + 3.818e-7 * exp(-1141. / T) + 2.064e-7 * exp(-3210. / T)
                  + 2.71e-5 * exp(-172700. / T) + 0.003802 * exp(-452300. / T) + 0.003086 * exp(-910500. / T));
    }
    inline double alphaDRMgIII(double T)
    {
        return pow(T, -1.5)
               * (6.269e-6 * exp(-410400. / T) + 0.0009181 * exp(-576600. / T) + 0.0003082 * exp(-731000. / T));
    }
    inline double alphaDRMgII(double T)
    {
        return pow(T, -1.5)
               * (3.871e-8 * exp(-8415. / T) + 4.732e-7 * exp(-16820. / T) + 0.001599 * exp(-50000. / T)
                  + 2.628e-5 * exp(-275900. / T));
    }

    // -- collisional ionization (Voronov 1997) --

    inline double betaMgI(double T)
    {
        double U = VernerCrossSections::MgI_eV * eV_to_K / T;
        return 6.21e-7 * pow(U, 0.39) * exp(-U) / (0.592 + U);
    }
    inline double betaMgII(double T)
    {
        double U = VernerCrossSections::MgII_eV * eV_to_K / T;
        return 1.92e-8 * pow(U, 0.85) * exp(-U) / (0.0027 + U);
    }
    inline double betaMgIII(double T)
    {
        double U = VernerCrossSections::MgIII_eV * eV_to_K / T;
        return 5.56e-9 * (1. + sqrt(U)) * pow(U, 0.3) * exp(-U) / (0.107 + U);
    }
    inline double betaMgIV(double T)
    {
        double U = VernerCrossSections::MgIV_eV * eV_to_K / T;
        return 4.35e-9 * (1. + sqrt(U)) * pow(U, 0.31) * exp(-U) / (0.159 + U);
    }
    inline double betaMgV(double T)
    {
        double U = VernerCrossSections::MgV_eV * eV_to_K / T;
        return 7.1e-9 * pow(U, 0.25) * exp(-U) / (0.658 + U);
    }
    inline double betaMgVI(double T)
    {
        double U = VernerCrossSections::MgVI_eV * eV_to_K / T;
        return 1.7e-9 * (1. + sqrt(U)) * pow(U, 0.28) * exp(-U) / (0.242 + U);
    }
    inline double betaMgVII(double T)
    {
        double U = VernerCrossSections::MgVII_eV * eV_to_K / T;
        return 1.22e-9 * (1. + sqrt(U)) * pow(U, 0.23) * exp(-U) / (0.343 + U);
    }
    inline double betaMgVIII(double T)
    {
        double U = VernerCrossSections::MgVIII_eV * eV_to_K / T;
        return 2.2e-9 * pow(U, 0.22) * exp(-U) / (0.897 + U);
    }
    inline double betaMgIX(double T)
    {
        double U = VernerCrossSections::MgIX_eV * eV_to_K / T;
        return 4.86e-10 * (1. + sqrt(U)) * pow(U, 0.14) * exp(-U) / (0.751 + U);
    }
    inline double betaMgX(double T)
    {
        double U = VernerCrossSections::MgX_eV * eV_to_K / T;
        return 2.35e-10 * (1. + sqrt(U)) * pow(U, 0.1) * exp(-U) / (1.03 + U);
    }

    // -- charge exchange --

    inline double xiHIMgV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 3e4) ? 3. : 1e-4 * T);
        return 6.36e-09 * pow(T4, 0.55) * (1. + 3.86 * exp(-5.19 * T4));
    }
    inline double xiHIMgIV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 3e4) ? 3. : 1e-4 * T);
        return 6.49e-09 * pow(T4, 0.53) * (1. + 2.82 * exp(-7.63 * T4));
    }
    inline double xiHIMgIII(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 3e4) ? 3. : 1e-4 * T);
        return 8.58e-14 * pow(T4, 0.00249) * (1. + 0.0293 * exp(-4.33 * T4));
    }
    inline double chiHIIMgII(double T)
    {
        double T4 = (T < 1e4) ? 1. : ((T > 3e5) ? 30. : 1e-4 * T);
        return 7.6e-14 * (1. - 1.97 * exp(-4.32 * T4)) * exp(-1.67 / T4);
    }
    inline double chiHIIMgI(double T)
    {
        double T4 = (T < 5e3) ? 0.5 : ((T > 3e4) ? 3. : 1e-4 * T);
        return 9.76e-12 * pow(T4, 3.14) * (1. + 55.54 * exp(-1.12 * T4));
    }

    // ============================================================
    //  Silicon  (SiI-SiXII: 12 ionization stages, rates up to SiXIII)
    // ============================================================

    // -- radiative recombination (ion+1 -> ion) --

    inline double alphaRRSiXIII(double T)
    {
        double s0 = sqrt(T / 6494.), s1 = sqrt(T / 5.693e7);
        return 2.017e-10 / (s0 * pow(1. + s0, 0.4412) * pow(1. + s1, 1.5588));
    }
    inline double alphaRRSiXII(double T)
    {
        double s0 = sqrt(T / 1088.), s1 = sqrt(T / 2.896e7);
        return 4.633e-10 / (s0 * pow(1. + s0, 0.3398) * pow(1. + s1, 1.6602));
    }
    inline double alphaRRSiXI(double T)
    {
        double s0 = sqrt(T / 69.06), s1 = sqrt(T / 1.644e7);
        return 1.851e-9 / (s0 * pow(1. + s0, 0.2616) * pow(1. + s1, 1.7384));
    }
    inline double alphaRRSiX(double T)
    {
        double s0 = sqrt(T / 55.49), s1 = sqrt(T / 1.716e7);
        return 1.688e-9 / (s0 * pow(1. + s0, 0.261) * pow(1. + s1, 1.739));
    }
    inline double alphaRRSiIX(double T)
    {
        double s0 = sqrt(T / 25.23), s1 = sqrt(T / 1.842e7);
        return 2.1e-9 / (s0 * pow(1. + s0, 0.2599) * pow(1. + s1, 1.7401));
    }
    inline double alphaRRSiVIII(double T)
    {
        double s0 = sqrt(T / 88.6), s1 = sqrt(T / 1.997e7);
        double Bp = 0.7072 + 0.0185 * exp(-694900. / T);
        return 7.532e-10 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRSiVII(double T)
    {
        double s0 = sqrt(T / 114.3), s1 = sqrt(T / 2.377e7);
        double Bp = 0.6753 + 0.0356 * exp(-859500. / T);
        return 4.615e-10 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRSiVI(double T)
    {
        double s0 = sqrt(T / 164.9), s1 = sqrt(T / 3.231e7);
        double Bp = 0.6113 + 0.0636 * exp(-983700. / T);
        return 2.468e-10 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRSiV(double T)
    {
        double s0 = sqrt(T / 1009.), s1 = sqrt(T / 8.514e7);
        double Bp = 0.3678 + 0.1646 * exp(-1.084e6 / T);
        return 5.134e-11 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRSiIV(double T)
    {
        double s0 = sqrt(T / 216.6), s1 = sqrt(T / 4.491e7);
        double Bp = 0.4931 + 0.1667 * exp(-904600. / T);
        return 6.739e-11 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRSiIII(double T)
    {
        double s0 = sqrt(T / 7.712), s1 = sqrt(T / 2.951e7);
        double Bp = 0.6287 + 0.1523 * exp(-480400. / T);
        return 1.964e-10 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRSiII(double T)
    {
        double s0 = sqrt(T / 15.9), s1 = sqrt(T / 4.237e7);
        double Bp = 0.627 + 0.2333 * exp(-58280. / T);
        return 3.262e-11 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }

    // -- dielectronic recombination --

    inline double alphaDRSiXIII(double T)
    {
        return pow(T, -1.5)
               * (0.05318 * exp(-1.552e7 / T) + 0.1874 * exp(-1.969e7 / T) + 0.01227 * exp(-2.532e7 / T)
                  + 0.0007173 * exp(-2.696e8 / T));
    }
    inline double alphaDRSiXII(double T)
    {
        return pow(T, -1.5)
               * (0.001205 * exp(-48240. / T) + 0.01309 * exp(-213700. / T) + 0.005333 * exp(-549200. / T)
                  + 0.02858 * exp(-2.46e6 / T) + 0.03195 * exp(-3.647e6 / T) + 0.1433 * exp(-1.889e7 / T));
    }
    inline double alphaDRSiXI(double T)
    {
        return pow(T, -1.5)
               * (0.0002214 * exp(-3586. / T) + 0.00126 * exp(-13240. / T) + 0.003832 * exp(-56360. / T)
                  + 0.01158 * exp(-249800. / T) + 0.02871 * exp(-536300. / T) + 0.05456 * exp(-2.809e6 / T)
                  + 0.05035 * exp(-1.854e7 / T));
    }
    inline double alphaDRSiX(double T)
    {
        return pow(T, -1.5)
               * (0.0001246 * exp(-1167. / T) + 0.0006649 * exp(-9088. / T) + 0.002912 * exp(-53320. / T)
                  + 0.02912 * exp(-421900. / T) + 0.03049 * exp(-1.571e6 / T) + 0.1056 * exp(-2.721e6 / T));
    }
    inline double alphaDRSiIX(double T)
    {
        return pow(T, -1.5)
               * (0.0005845 * exp(-985.6 / T) + 0.00066 * exp(-6577. / T) + 0.000718 * exp(-42810. / T)
                  + 0.01714 * exp(-394400. / T) + 0.02958 * exp(-1.296e6 / T) + 0.1075 * exp(-2.441e6 / T));
    }
    inline double alphaDRSiVIII(double T)
    {
        return pow(T, -1.5)
               * (5.272e-5 * exp(-2829. / T) + 0.0002282 * exp(-26170. / T) + 0.001345 * exp(-137400. / T)
                  + 0.02246 * exp(-452000. / T) + 0.09606 * exp(-2.072e6 / T) + 0.0008366 * exp(-7.808e6 / T));
    }
    inline double alphaDRSiVII(double T)
    {
        return pow(T, -1.5)
               * (2.086e-6 * exp(-370.8 / T) + 9.423e-6 * exp(-3870. / T) + 3.423e-5 * exp(-22260. / T)
                  + 0.000395 * exp(-131800. / T) + 0.01535 * exp(-528500. / T) + 0.04986 * exp(-1.742e6 / T)
                  + 0.0004067 * exp(-8.392e6 / T));
    }
    inline double alphaDRSiVI(double T)
    {
        return pow(T, -1.5)
               * (7.163e-7 * exp(-562.5 / T) + 2.656e-6 * exp(-2952. / T) + 1.119e-6 * exp(-9682. / T)
                  + 4.796e-5 * exp(-147300. / T) + 0.004052 * exp(-506400. / T) + 0.006101 * exp(-804700. / T)
                  + 0.02366 * exp(-1.623e6 / T));
    }
    inline double alphaDRSiV(double T)
    {
        return pow(T, -1.5)
               * (0.0001422 * exp(-768500. / T) + 0.009474 * exp(-1.208e6 / T) + 0.00165 * exp(-1.839e6 / T));
    }
    inline double alphaDRSiIV(double T)
    {
        return pow(T, -1.5)
               * (3.819e-6 * exp(-3802. / T) + 2.421e-5 * exp(-12800. / T) + 0.0002283 * exp(-59530. / T)
                  + 0.008604 * exp(-102600. / T) + 0.002617 * exp(-1.154e6 / T));
    }
    inline double alphaDRSiIII(double T)
    {
        return pow(T, -1.5)
               * (2.93e-6 * exp(-116.2 / T) + 2.803e-6 * exp(-5721. / T) + 9.023e-5 * exp(-34770. / T)
                  + 0.006909 * exp(-117600. / T) + 2.582e-5 * exp(-3.505e6 / T));
    }
    inline double alphaDRSiII(double T)
    {
        return pow(T, -1.5)
               * (3.408e-8 * exp(-24.31 / T) + 1.913e-7 * exp(-129.3 / T) + 1.679e-7 * exp(-427.2 / T)
                  + 7.523e-7 * exp(-3729. / T) + 8.386e-5 * exp(-55140. / T) + 0.004083 * exp(-129500. / T));
    }

    // -- collisional ionization (Voronov 1997) --

    inline double betaSiI(double T)
    {
        double U = VernerCrossSections::SiI_eV * eV_to_K / T;
        return 1.88e-7 * (1. + sqrt(U)) * pow(U, 0.25) * exp(-U) / (0.376 + U);
    }
    inline double betaSiII(double T)
    {
        double U = VernerCrossSections::SiII_eV * eV_to_K / T;
        return 6.43e-8 * (1. + sqrt(U)) * pow(U, 0.2) * exp(-U) / (0.632 + U);
    }
    inline double betaSiIII(double T)
    {
        double U = VernerCrossSections::SiIII_eV * eV_to_K / T;
        return 2.01e-8 * (1. + sqrt(U)) * pow(U, 0.22) * exp(-U) / (0.473 + U);
    }
    inline double betaSiIV(double T)
    {
        double U = VernerCrossSections::SiIV_eV * eV_to_K / T;
        return 4.94e-9 * (1. + sqrt(U)) * pow(U, 0.23) * exp(-U) / (0.172 + U);
    }
    inline double betaSiV(double T)
    {
        double U = VernerCrossSections::SiV_eV * eV_to_K / T;
        return 1.76e-9 * (1. + sqrt(U)) * pow(U, 0.31) * exp(-U) / (0.102 + U);
    }
    inline double betaSiVI(double T)
    {
        double U = VernerCrossSections::SiVI_eV * eV_to_K / T;
        return 1.74e-9 * (1. + sqrt(U)) * pow(U, 0.29) * exp(-U) / (0.18 + U);
    }
    inline double betaSiVII(double T)
    {
        double U = VernerCrossSections::SiVII_eV * eV_to_K / T;
        return 1.23e-9 * (1. + sqrt(U)) * pow(U, 0.07) * exp(-U) / (0.518 + U);
    }
    inline double betaSiVIII(double T)
    {
        double U = VernerCrossSections::SiVIII_eV * eV_to_K / T;
        return 8.27e-10 * (1. + sqrt(U)) * pow(U, 0.28) * exp(-U) / (0.239 + U);
    }
    inline double betaSiIX(double T)
    {
        double U = VernerCrossSections::SiIX_eV * eV_to_K / T;
        return 6.01e-10 * (1. + sqrt(U)) * pow(U, 0.25) * exp(-U) / (0.305 + U);
    }
    inline double betaSiX(double T)
    {
        double U = VernerCrossSections::SiX_eV * eV_to_K / T;
        return 4.65e-10 * (1. + sqrt(U)) * pow(U, 0.04) * exp(-U) / (0.666 + U);
    }
    inline double betaSiXI(double T)
    {
        double U = VernerCrossSections::SiXI_eV * eV_to_K / T;
        return 2.63e-10 * (1. + sqrt(U)) * pow(U, 0.16) * exp(-U) / (0.666 + U);
    }
    inline double betaSiXII(double T)
    {
        double U = VernerCrossSections::SiXII_eV * eV_to_K / T;
        return 1.18e-10 * (1. + sqrt(U)) * pow(U, 0.16) * exp(-U) / (0.734 + U);
    }

    // -- charge exchange --

    inline double xiHISiV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 5e4) ? 5. : 1e-4 * T);
        return 7.58e-9 * pow(T4, 0.37) * (1. + 1.06 * exp(-4.09 * T4));
    }
    inline double xiHISiIV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 3e4) ? 3. : 1e-4 * T);
        return 4.9e-10 * pow(T4, -0.0874) * (1. - 0.36 * exp(-0.79 * T4));
    }
    inline double xiHISiIII(double T)
    {
        double T4 = (T < 500.) ? 0.05 : ((T > 1e5) ? 10. : 1e-4 * T);
        return 6.77e-9 * pow(T4, 0.0736) * (1. - 0.43 * exp(-0.11 * T4));
    }
    inline double chiHIISiII(double T)
    {
        double T4 = (T < 2e3) ? 0.2 : ((T > 1e5) ? 10. : 1e-4 * T);
        return 2.26e-9 * pow(T4, 0.0736) * (1. - 0.43 * exp(-0.11 * T4)) * exp(-3.031 / T4);
    }
    inline double chiHIISiI(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 2e5) ? 20. : 1e-4 * T);
        return 9.2e-10 * pow(T4, 1.15) * (1. + 0.8 * exp(-0.24 * T4));
    }

    // ============================================================
    //  Sulfur  (SI-SIV: first 4 ionization stages, matching CHIANTI)
    // ============================================================

    inline double alphaRRSV(double T)
    {
        double s0 = sqrt(T / 62.38), s1 = sqrt(T / 2.803e7);
        double Bp = 0.6343 + 0.0773 * exp(-1.059e6 / T);
        return 2.615e-10 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRSIV(double T)
    {
        double s0 = sqrt(T / 16.78), s1 = sqrt(T / 2.05e7);
        double Bp = 0.6947 + 0.0795 * exp(-68680. / T);
        return 3.043e-10 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRSIII(double T)
    {
        double s0 = sqrt(T / 329.4), s1 = sqrt(T / 2.166e7);
        double Bp = 0.4642 + 0.3351 * exp(-763000. / T);
        return 2.478e-11 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRSII(double T)
    {
        return 4.1e-13 * pow(1e-4 * T, -0.63);
    }

    inline double alphaDRSV(double T)
    {
        return pow(T, -1.5)
               * (9.571e-6 * exp(-1180. / T) + 6.268e-5 * exp(-6443. / T) + 0.0003807 * exp(-22640. / T)
                  + 0.01874 * exp(-153000. / T) + 0.005526 * exp(-356400. / T));
    }
    inline double alphaDRSIV(double T)
    {
        return pow(T, -1.5)
               * (5.817e-7 * exp(-362.8 / T) + 1.391e-6 * exp(-1058. / T) + 1.123e-5 * exp(-7160. / T)
                  + 0.0001521 * exp(-32600. / T) + 0.001875 * exp(-123500. / T) + 0.02097 * exp(-207000. / T));
    }
    inline double alphaDRSIII(double T)
    {
        return pow(T, -1.5)
               * (3.04e-7 * exp(-50.16 / T) + 4.393e-7 * exp(-326.6 / T) + 1.609e-6 * exp(-3102. / T)
                  + 4.98e-6 * exp(-12100. / T) + 3.457e-5 * exp(-49690. / T) + 0.008617 * exp(-201000. / T)
                  + 0.0009284 * exp(-257500. / T));
    }
    inline double alphaDRSII(double T)
    {
        return 0.0017126 * pow(T, -1.5) * exp(-173490. / T);
    }

    inline double betaSI(double T)
    {
        double U = VernerCrossSections::SI_eV * eV_to_K / T;
        return 5.49e-8 * (1. + sqrt(U)) * pow(U, 0.25) * exp(-U) / (0.1 + U);
    }
    inline double betaSII(double T)
    {
        double U = VernerCrossSections::SII_eV * eV_to_K / T;
        return 6.81e-8 * (1. + sqrt(U)) * pow(U, 0.21) * exp(-U) / (0.693 + U);
    }
    inline double betaSIII(double T)
    {
        double U = VernerCrossSections::SIII_eV * eV_to_K / T;
        return 2.14e-8 * (1. + sqrt(U)) * pow(U, 0.24) * exp(-U) / (0.353 + U);
    }
    inline double betaSIV(double T)
    {
        double U = VernerCrossSections::SIV_eV * eV_to_K / T;
        return 1.66e-8 * (1. + sqrt(U)) * pow(U, 0.14) * exp(-U) / (1.03 + U);
    }

    inline double xiHISIV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 3e4) ? 3. : 1e-4 * T);
        return 2.29e-9 * pow(T4, 0.0402) * (1. + 1.59 * exp(-6.06 * T4));
    }
    inline double xiHISIII(double /*T*/)
    {
        return 1e-14;
    }
    inline double xiHISII(double T)
    {
        double T4 = (T < 3e3) ? 0.3 : ((T > 1e4) ? 1. : 1e-4 * T);
        return 3.82e-16 * pow(T4, 11.1) * (1. + 25700. * exp(-8.22 * T4));
    }
    inline double chiHIISI(double /*T*/)
    {
        return 1e-14;
    }

    // ============================================================
    //  Iron  (FeI-FeXVI: 16 ionization stages -> FeXVII is bare nucleus)
    // ============================================================

    inline double alphaRRFeXVII(double T)
    {
        double s0 = sqrt(T / 17510.), s1 = sqrt(T / 1.579e8);
        return 2.034e-10 / (s0 * pow(1. + s0, 0.5452) * pow(1. + s1, 1.4548));
    }
    inline double alphaRRFeXVI(double T)
    {
        double s0 = sqrt(T / 6295.), s1 = sqrt(T / 9.035e7);
        return 3.133e-10 / (s0 * pow(1. + s0, 0.4493) * pow(1. + s1, 1.5507));
    }
    inline double alphaRRFeXV(double T)
    {
        double s0 = sqrt(T / 1881.), s1 = sqrt(T / 5.429e7);
        return 5.398e-10 / (s0 * pow(1. + s0, 0.3705) * pow(1. + s1, 1.6295));
    }
    inline double alphaRRFeXIV(double T)
    {
        double s0 = sqrt(T / 1531.), s1 = sqrt(T / 4.632e7);
        double Bp = 0.6338 + 0.0278 * exp(-90700. / T);
        return 5.37e-10 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRFeXIII(double T)
    {
        double s0 = sqrt(T / 115.8), s1 = sqrt(T / 4.4e7);
        return 1.984e-9 / (s0 * pow(1. + s0, 0.2899) * pow(1. + s1, 1.7101));
    }
    inline double alphaRRFeXII(double T)
    {
        double s0 = sqrt(T / 353.1), s1 = sqrt(T / 3.554e7);
        double Bp = 0.7156 + 0.0132 * exp(-295100. / T);
        return 8.303e-10 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRFeXI(double T)
    {
        double s0 = sqrt(T / 163.9), s1 = sqrt(T / 2.924e7);
        double Bp = 0.737 + 0.0224 * exp(-429100. / T);
        return 1.052e-9 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRFeX(double T)
    {
        double s0 = sqrt(T / 72.42), s1 = sqrt(T / 2.453e7);
        double Bp = 0.7495 + 0.0404 * exp(-419900. / T);
        return 1.338e-9 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRFeIX(double T)
    {
        double s0 = sqrt(T / 1891.), s1 = sqrt(T / 5.181e7);
        double Bp = 0.4865 + 0.0575 * exp(-27340. / T);
        return 3.341e-10 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRFeVIII(double T)
    {
        double s0 = sqrt(T / 37.53), s1 = sqrt(T / 2.418e7);
        double Bp = 0.7407 + 0.072 * exp(-328000. / T);
        return 1.121e-9 / (s0 * pow(1. + s0, 1. - Bp) * pow(1. + s1, 1. + Bp));
    }
    inline double alphaRRFeVII(double T)
    {
        return 2.62e-11 * pow(1e-4 * T, -0.728);
    }
    inline double alphaRRFeVI(double T)
    {
        return 1.51e-11 * pow(1e-4 * T, -0.699);
    }
    inline double alphaRRFeV(double T)
    {
        return 7.8e-12 * pow(1e-4 * T, -0.682);
    }
    inline double alphaRRFeIV(double T)
    {
        return 3.32e-12 * pow(1e-4 * T, -0.746);
    }
    inline double alphaRRFeIII(double T)
    {
        return 1.02e-12 * pow(1e-4 * T, -0.843);
    }
    inline double alphaRRFeII(double T)
    {
        return 1.42e-13 * pow(1e-4 * T, -0.891);
    }

    inline double alphaDRFeXVI(double T)
    {
        return pow(T, -1.5)
               * (0.0007676 * exp(-29350. / T) + 0.005587 * exp(-81580. / T) + 0.1152 * exp(-359100. / T)
                  + 0.04929 * exp(-1.735e6 / T) + 0.7274 * exp(-7.545e6 / T) + 0.007347 * exp(-4.634e7 / T));
    }
    inline double alphaDRFeXV(double T)
    {
        return pow(T, -1.5)
               * (0.0005636 * exp(-3628. / T) + 0.00786 * exp(-24890. / T) + 0.05063 * exp(-140500. / T)
                  + 0.1753 * exp(-513300. / T) + 0.1209 * exp(-5.018e6 / T) + 0.1934 * exp(-8.689e6 / T));
    }
    inline double alphaDRFeXIV(double T)
    {
        return pow(T, -1.5)
               * (0.002305 * exp(-4747. / T) + 0.01072 * exp(-18770. / T) + 0.03512 * exp(-119000. / T)
                  + 0.2105 * exp(-509000. / T) + 0.03622 * exp(-1.595e6 / T));
    }
    inline double alphaDRFeXIII(double T)
    {
        return pow(T, -1.5)
               * (0.004469 * exp(-2462. / T) + 0.008538 * exp(-12610. / T) + 0.01741 * exp(-93300. / T)
                  + 0.163 * exp(-488700. / T) + 0.0868 * exp(-1.312e6 / T));
    }
    inline double alphaDRFeXII(double T)
    {
        return pow(T, -1.5)
               * (0.001074 * exp(-1387. / T) + 0.00608 * exp(-10480. / T) + 0.01887 * exp(-39550. / T)
                  + 0.0254 * exp(-146100. / T) + 0.0758 * exp(-401000. / T) + 0.2773 * exp(-720800. / T));
    }
    inline double alphaDRFeXI(double T)
    {
        return pow(T, -1.5)
               * (6.487e-5 * exp(-110.1 / T) + 8.793e-5 * exp(-565.4 / T) + 0.0004939 * exp(-1842. / T)
                  + 0.003787 * exp(-7134. / T) + 0.008878 * exp(-30850. / T) + 0.05325 * exp(-187800. / T)
                  + 0.2104 * exp(-670600. / T));
    }
    inline double alphaDRFeX(double T)
    {
        return pow(T, -1.5)
               * (6.485e-5 * exp(-39.94 / T) + 6.36e-5 * exp(-562.1 / T) + 0.000372 * exp(-1992. / T)
                  + 0.001607 * exp(-8325. / T) + 0.003516 * exp(-27570. / T) + 0.007326 * exp(-74090. / T)
                  + 0.0256 * exp(-155200. / T) + 0.1005 * exp(-438800. / T) + 0.1942 * exp(-735500. / T));
    }
    inline double alphaDRFeIX(double T)
    {
        return pow(T, -1.5)
               * (4.777e-7 * exp(-9.034 / T) + 1.231e-6 * exp(-112.8 / T) + 5.055e-5 * exp(-662.4 / T)
                  + 0.0003413 * exp(-1143. / T) + 0.001625 * exp(-3926. / T) + 0.003873 * exp(-13000. / T)
                  + 0.006438 * exp(-46840. / T) + 0.0697 * exp(-267000. / T) + 0.2925 * exp(-735800. / T));
    }
    inline double alphaDRFeVIII(double T)
    {
        return pow(T, -1.5)
               * (5.978e-7 * exp(-8.385 / T) + 8.939e-7 * exp(-99.22 / T) + 1.64e-5 * exp(-523.4 / T)
                  + 9.598e-5 * exp(-1579. / T) + 0.0001105 * exp(-4489. / T) + 0.0007299 * exp(-21020. / T)
                  + 0.003858 * exp(-97780. / T) + 0.02476 * exp(-335300. / T) + 0.1789 * exp(-808100. / T));
    }
    inline double alphaDRFeVII(double T)
    {
        return pow(T, -1.5) * (0.092054 * exp(-528000. / T) + 0.041024 * exp(-4.1776e6 / T));
    }
    inline double alphaDRFeVI(double T)
    {
        return pow(T, -1.5) * (0.080047 * exp(-628960. / T) + 0.024014 * exp(-1.1605e6 / T));
    }
    inline double alphaDRFeV(double T)
    {
        return pow(T, -1.5) * (0.038023 * exp(-432850. / T) + 0.01601 * exp(-782140. / T));
    }
    inline double alphaDRFeIV(double T)
    {
        return pow(T, -1.5) * (0.015009 * exp(-331890. / T) + 0.0047027 * exp(-604590. / T));
    }
    inline double alphaDRFeIII(double T)
    {
        return pow(T, -1.5) * (0.0023013 * exp(-193800. / T) + 0.0027016 * exp(-364380. / T));
    }
    inline double alphaDRFeII(double T)
    {
        return pow(T, -1.5) * (0.00022013 * exp(-59415. / T) + 0.00010006 * exp(-149700. / T));
    }

    inline double betaFeI(double T)
    {
        double U = VernerCrossSections::FeI_eV * eV_to_K / T;
        return 2.52e-7 * pow(U, 0.25) * exp(-U) / (0.701 + U);
    }
    inline double betaFeII(double T)
    {
        double U = VernerCrossSections::FeII_eV * eV_to_K / T;
        return 2.21e-8 * (1. + sqrt(U)) * pow(U, 0.45) * exp(-U) / (0.033 + U);
    }
    inline double betaFeIII(double T)
    {
        double U = VernerCrossSections::FeIII_eV * eV_to_K / T;
        return 4.1e-8 * pow(U, 0.17) * exp(-U) / (0.366 + U);
    }
    inline double betaFeIV(double T)
    {
        double U = VernerCrossSections::FeIV_eV * eV_to_K / T;
        return 3.53e-8 * pow(U, 0.39) * exp(-U) / (0.243 + U);
    }
    inline double betaFeV(double T)
    {
        double U = VernerCrossSections::FeV_eV * eV_to_K / T;
        return 1.04e-8 * (1. + sqrt(U)) * pow(U, 0.17) * exp(-U) / (0.285 + U);
    }
    inline double betaFeVI(double T)
    {
        double U = VernerCrossSections::FeVI_eV * eV_to_K / T;
        return 1.23e-8 * (1. + sqrt(U)) * pow(U, 0.21) * exp(-U) / (0.411 + U);
    }
    inline double betaFeVII(double T)
    {
        double U = VernerCrossSections::FeVII_eV * eV_to_K / T;
        return 9.47e-9 * (1. + sqrt(U)) * pow(U, 0.21) * exp(-U) / (0.458 + U);
    }
    inline double betaFeVIII(double T)
    {
        double U = VernerCrossSections::FeVIII_eV * eV_to_K / T;
        return 4.71e-9 * (1. + sqrt(U)) * pow(U, 0.28) * exp(-U) / (0.28 + U);
    }
    inline double betaFeIX(double T)
    {
        double U = VernerCrossSections::FeIX_eV * eV_to_K / T;
        return 3.02e-9 * (1. + sqrt(U)) * pow(U, 0.15) * exp(-U) / (0.697 + U);
    }
    inline double betaFeX(double T)
    {
        double U = VernerCrossSections::FeX_eV * eV_to_K / T;
        return 2.34e-9 * (1. + sqrt(U)) * pow(U, 0.14) * exp(-U) / (0.764 + U);
    }
    inline double betaFeXI(double T)
    {
        double U = VernerCrossSections::FeXI_eV * eV_to_K / T;
        return 1.76e-9 * (1. + sqrt(U)) * pow(U, 0.14) * exp(-U) / (0.805 + U);
    }
    inline double betaFeXII(double T)
    {
        double U = VernerCrossSections::FeXII_eV * eV_to_K / T;
        return 1.14e-9 * (1. + sqrt(U)) * pow(U, 0.15) * exp(-U) / (0.773 + U);
    }
    inline double betaFeXIII(double T)
    {
        double U = VernerCrossSections::FeXIII_eV * eV_to_K / T;
        return 8.66e-10 * (1. + sqrt(U)) * pow(U, 0.14) * exp(-U) / (0.805 + U);
    }
    inline double betaFeXIV(double T)
    {
        double U = VernerCrossSections::FeXIV_eV * eV_to_K / T;
        return 6.61e-10 * (1. + sqrt(U)) * pow(U, 0.14) * exp(-U) / (0.762 + U);
    }
    inline double betaFeXV(double T)
    {
        double U = VernerCrossSections::FeXV_eV * eV_to_K / T;
        return 4.41e-10 * (1. + sqrt(U)) * pow(U, 0.16) * exp(-U) / (0.698 + U);
    }
    inline double betaFeXVI(double T)
    {
        double U = VernerCrossSections::FeXVI_eV * eV_to_K / T;
        return 1.18e-10 * (1. + sqrt(U)) * pow(U, 0.15) * exp(-U) / (0.211 + U);
    }

    inline double xiHIFeV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 3e4) ? 3. : 1e-4 * T);
        return 1.46e-8 * pow(T4, 0.0357) * (1. - 0.92 * exp(-0.37 * T4));
    }
    inline double xiHIFeIV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 1e5) ? 10. : 1e-4 * T);
        return 3.42e-9 * pow(T4, 0.51) * (1. - 2.06 * exp(-8.99 * T4));
    }
    inline double xiHIFeIII(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 1e5) ? 10. : 1e-4 * T);
        return 1.26e-9 * pow(T4, 0.0772) * (1. - 0.41 * exp(-7.31 * T4));
    }
    inline double chiHIIFeII(double T)
    {
        double T4 = (T < 1e4) ? 1. : ((T > 1e5) ? 10. : 1e-4 * T);
        return 2.1e-9 * pow(T4, 0.0772) * (1. - 0.41 * exp(-7.31 * T4)) * exp(-3.005 / T4);
    }
    inline double chiHIIFeI(double /*T*/)
    {
        return 5.4e-9;
    }

    // ============================================================
    //  Dispatch functions by linear ion index
    //  (same ordering as VernerCrossSections)
    // ============================================================

    /** Total recombination rate alpha_RR + alpha_DR for the recombination from ion stage m+1 -> m.
        The ion index refers to the LOWER state (the state being recombined INTO).
        For H and He, returns Case B recombination. For metals, returns RR + DR. */
    inline double alphaTotal(int ion, double T)
    {
        switch (ion)
        {
            // H
            case 0: return alphaBHII(T);
            // He
            case 1: return alphaBHeII(T) + alphaDRHeII(T);
            case 2: return alphaBHeIII(T);
            // C  (ion 3=CI means CVII->CVI recomb uses alphaRRCVII, etc.)
            // Note: ion index = lower state. Recomb from (ion+1) -> ion.
            // For CI (index 3): recomb CII -> CI uses alphaRRCII + alphaDRCII
            case 3: return alphaRRCII(T) + alphaDRCII(T);
            case 4: return alphaRRCIII(T) + alphaDRCIII(T);
            case 5: return alphaRRCIV(T) + alphaDRCIV(T);
            case 6: return alphaRRCV(T) + alphaDRCV(T);
            case 7: return alphaRRCVI(T) + alphaDRCVI(T);
            case 8:
                return alphaRRCVII(T);  // bare C nucleus, no DR
            // N
            case 9: return alphaRRNII(T) + alphaDRNII(T);
            case 10: return alphaRRNIII(T) + alphaDRNIII(T);
            case 11: return alphaRRNIV(T) + alphaDRNIV(T);
            case 12: return alphaRRNV(T) + alphaDRNV(T);
            case 13: return alphaRRNVI(T) + alphaDRNVI(T);
            case 14: return alphaRRNVII(T) + alphaDRNVII(T);
            case 15: return alphaRRNVIII(T);
            // O
            case 16: return alphaRROII(T) + alphaDROII(T);
            case 17: return alphaRROIII(T) + alphaDROIII(T);
            case 18: return alphaRROIV(T) + alphaDROIV(T);
            case 19: return alphaRROV(T) + alphaDROV(T);
            case 20: return alphaRROVI(T) + alphaDROVI(T);
            case 21: return alphaRROVII(T) + alphaDROVII(T);
            case 22: return alphaRROVIII(T) + alphaDROVIII(T);
            case 23: return alphaRROIX(T);
            // Ne
            case 24: return alphaRRNeII(T) + alphaDRNeII(T);
            case 25: return alphaRRNeIII(T) + alphaDRNeIII(T);
            case 26: return alphaRRNeIV(T) + alphaDRNeIV(T);
            case 27: return alphaRRNeV(T) + alphaDRNeV(T);
            case 28: return alphaRRNeVI(T) + alphaDRNeVI(T);
            case 29: return alphaRRNeVII(T) + alphaDRNeVII(T);
            case 30: return alphaRRNeVIII(T) + alphaDRNeVIII(T);
            case 31: return alphaRRNeIX(T);
            // Mg
            case 32: return alphaRRMgII(T) + alphaDRMgII(T);
            case 33: return alphaRRMgIII(T) + alphaDRMgIII(T);
            case 34: return alphaRRMgIV(T) + alphaDRMgIV(T);
            case 35: return alphaRRMgV(T) + alphaDRMgV(T);
            case 36: return alphaRRMgVI(T) + alphaDRMgVI(T);
            case 37: return alphaRRMgVII(T) + alphaDRMgVII(T);
            case 38: return alphaRRMgVIII(T) + alphaDRMgVIII(T);
            case 39: return alphaRRMgIX(T) + alphaDRMgIX(T);
            case 40: return alphaRRMgX(T) + alphaDRMgX(T);
            case 41: return alphaRRMgXI(T) + alphaDRMgXI(T);
            // Si
            case 42: return alphaRRSiII(T) + alphaDRSiII(T);
            case 43: return alphaRRSiIII(T) + alphaDRSiIII(T);
            case 44: return alphaRRSiIV(T) + alphaDRSiIV(T);
            case 45: return alphaRRSiV(T) + alphaDRSiV(T);
            case 46: return alphaRRSiVI(T) + alphaDRSiVI(T);
            case 47: return alphaRRSiVII(T) + alphaDRSiVII(T);
            case 48: return alphaRRSiVIII(T) + alphaDRSiVIII(T);
            case 49: return alphaRRSiIX(T) + alphaDRSiIX(T);
            case 50: return alphaRRSiX(T) + alphaDRSiX(T);
            case 51: return alphaRRSiXI(T) + alphaDRSiXI(T);
            case 52: return alphaRRSiXII(T) + alphaDRSiXII(T);
            case 53: return alphaRRSiXIII(T) + alphaDRSiXIII(T);
            // S
            case 54: return alphaRRSII(T) + alphaDRSII(T);
            case 55: return alphaRRSIII(T) + alphaDRSIII(T);
            case 56: return alphaRRSIV(T) + alphaDRSIV(T);
            case 57: return alphaRRSV(T) + alphaDRSV(T);
            // Fe
            case 58: return alphaRRFeII(T) + alphaDRFeII(T);
            case 59: return alphaRRFeIII(T) + alphaDRFeIII(T);
            case 60: return alphaRRFeIV(T) + alphaDRFeIV(T);
            case 61: return alphaRRFeV(T) + alphaDRFeV(T);
            case 62: return alphaRRFeVI(T) + alphaDRFeVI(T);
            case 63: return alphaRRFeVII(T) + alphaDRFeVII(T);
            case 64: return alphaRRFeVIII(T) + alphaDRFeVIII(T);
            case 65: return alphaRRFeIX(T) + alphaDRFeIX(T);
            case 66: return alphaRRFeX(T) + alphaDRFeX(T);
            case 67: return alphaRRFeXI(T) + alphaDRFeXI(T);
            case 68: return alphaRRFeXII(T) + alphaDRFeXII(T);
            case 69: return alphaRRFeXIII(T) + alphaDRFeXIII(T);
            case 70: return alphaRRFeXIV(T) + alphaDRFeXIV(T);
            case 71: return alphaRRFeXV(T) + alphaDRFeXV(T);
            case 72: return alphaRRFeXVI(T) + alphaDRFeXVI(T);
            case 73: return alphaRRFeXVII(T);
            default: return 0.;
        }
    }

    /** Collisional ionization rate [cm^3/s] for ion stage m -> m+1.
        The ion index refers to the state being ionized FROM. */
    inline double beta(int ion, double T)
    {
        switch (ion)
        {
            case 0: return betaHI(T);
            case 1: return betaHeI(T);
            case 2: return betaHeII(T);
            case 3: return betaCI(T);
            case 4: return betaCII(T);
            case 5: return betaCIII(T);
            case 6: return betaCIV(T);
            case 7: return betaCV(T);
            case 8: return betaCVI(T);
            case 9: return betaNI(T);
            case 10: return betaNII(T);
            case 11: return betaNIII(T);
            case 12: return betaNIV(T);
            case 13: return betaNV(T);
            case 14: return betaNVI(T);
            case 15: return betaNVII(T);
            case 16: return betaOI(T);
            case 17: return betaOII(T);
            case 18: return betaOIII(T);
            case 19: return betaOIV(T);
            case 20: return betaOV(T);
            case 21: return betaOVI(T);
            case 22: return betaOVII(T);
            case 23: return betaOVIII(T);
            case 24: return betaNeI(T);
            case 25: return betaNeII(T);
            case 26: return betaNeIII(T);
            case 27: return betaNeIV(T);
            case 28: return betaNeV(T);
            case 29: return betaNeVI(T);
            case 30: return betaNeVII(T);
            case 31: return betaNeVIII(T);
            // Mg
            case 32: return betaMgI(T);
            case 33: return betaMgII(T);
            case 34: return betaMgIII(T);
            case 35: return betaMgIV(T);
            case 36: return betaMgV(T);
            case 37: return betaMgVI(T);
            case 38: return betaMgVII(T);
            case 39: return betaMgVIII(T);
            case 40: return betaMgIX(T);
            case 41: return betaMgX(T);
            // Si
            case 42: return betaSiI(T);
            case 43: return betaSiII(T);
            case 44: return betaSiIII(T);
            case 45: return betaSiIV(T);
            case 46: return betaSiV(T);
            case 47: return betaSiVI(T);
            case 48: return betaSiVII(T);
            case 49: return betaSiVIII(T);
            case 50: return betaSiIX(T);
            case 51: return betaSiX(T);
            case 52: return betaSiXI(T);
            case 53: return betaSiXII(T);
            // S
            case 54: return betaSI(T);
            case 55: return betaSII(T);
            case 56: return betaSIII(T);
            case 57: return betaSIV(T);
            // Fe
            case 58: return betaFeI(T);
            case 59: return betaFeII(T);
            case 60: return betaFeIII(T);
            case 61: return betaFeIV(T);
            case 62: return betaFeV(T);
            case 63: return betaFeVI(T);
            case 64: return betaFeVII(T);
            case 65: return betaFeVIII(T);
            case 66: return betaFeIX(T);
            case 67: return betaFeX(T);
            case 68: return betaFeXI(T);
            case 69: return betaFeXII(T);
            case 70: return betaFeXIII(T);
            case 71: return betaFeXIV(T);
            case 72: return betaFeXV(T);
            case 73: return betaFeXVI(T);
            default: return 0.;
        }
    }

    // ============================================================
    //  Charge exchange dispatch by linear ion index
    // ============================================================

    /** Hydrogen charge exchange recombination rate [cm^3/s] for a given ion (Verner index).
        Reaction: ion + HI -> ion_lower + HII.
        The ion index refers to the HIGHER stage (the one being recombined). */
    inline double xiHI(int ion, double T)
    {
        switch (ion)
        {
            // He
            case 1: return xiHIHeII(T);
            case 2: return xiHIHeIII(T);
            // C: CII-CV (ions 4-7 in Verner indexing)
            case 4: return xiHICII(T);
            case 5: return xiHICIII(T);
            case 6: return xiHICIV(T);
            case 7: return xiHICV(T);
            // N: NII-NV (ions 10-13)
            case 10: return xiHINII(T);
            case 11: return xiHINIII(T);
            case 12: return xiHINIV(T);
            case 13: return xiHINV(T);
            // O: OII-OV (ions 17-20)
            case 17: return xiHIOII(T);
            case 18: return xiHIOIII(T);
            case 19: return xiHIOIV(T);
            case 20: return xiHIOV(T);
            // Ne: NeIII-NeV (ions 26-28)
            case 26: return xiHINeIII(T);
            case 27: return xiHINeIV(T);
            case 28: return xiHINeV(T);
            // Mg: MgIII-MgV (ions 34-36)
            case 34: return xiHIMgIII(T);
            case 35: return xiHIMgIV(T);
            case 36: return xiHIMgV(T);
            // Si: SiIII-SiV (ions 44-46)
            case 44: return xiHISiIII(T);
            case 45: return xiHISiIV(T);
            case 46: return xiHISiV(T);
            // S: SII-SIV (ions 55-57)
            case 55: return xiHISII(T);
            case 56: return xiHISIII(T);
            case 57: return xiHISIV(T);
            // Fe: FeIII-FeV (ions 60-62)
            case 60: return xiHIFeIII(T);
            case 61: return xiHIFeIV(T);
            case 62: return xiHIFeV(T);
            default: return 0.;
        }
    }

    /** Hydrogen charge exchange ionization rate [cm^3/s] for a given ion (Verner index).
        Reaction: ion + HII -> ion_higher + HI.
        The ion index refers to the LOWER stage (the one being ionized). */
    inline double chiHII(int ion, double T)
    {
        switch (ion)
        {
            // C: CI (ion 3)
            case 3: return chiHIICI(T);
            // N: NI (ion 9)
            case 9: return chiHIINI(T);
            // O: OI (ion 16)
            case 16: return chiHIIOI(T);
            // Mg: MgI-MgII (ions 32-33)
            case 32: return chiHIIMgI(T);
            case 33: return chiHIIMgII(T);
            // Si: SiI-SiII (ions 42-43)
            case 42: return chiHIISiI(T);
            case 43: return chiHIISiII(T);
            // S: SI (ion 54)
            case 54: return chiHIISI(T);
            // Fe: FeI-FeII (ions 58-59)
            case 58: return chiHIIFeI(T);
            case 59: return chiHIIFeII(T);
            default: return 0.;
        }
    }

    // ============================================================
    //  Helium charge exchange rates (Kingdon & Ferland 1996)
    //  xi_HeI: ion+ + HeI -> ion + HeII (recombination)
    //  chi_HeII: ion + HeII -> ion+ + HeI (ionization)
    // ============================================================

    // -- C --
    inline double xiHeICIV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 1e7) ? 1000. : 1e-4 * T);
        return 1.12e-9 * pow(T4, 0.42) * (1. - 0.69 * exp(-0.34 * T4));
    }
    inline double xiHeICV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 1e7) ? 1000. : 1e-4 * T);
        if (T4 < 1.) return 3.12e-16 * pow(T4, -0.0737) * (1. + 35. * exp(2.4 * T4));
        if (T4 < 35.) return 1.49e-14 * pow(T4, 2.73) * (1. + 5.93 * exp(-0.0874 * T4));
        return 5.8e-11 * pow(T4, 0.73) * (1. - 0.86 * exp(-0.0096 * T4));
    }
    inline double chiHeIICI(double T)
    {
        return 6.3e-15 * pow(T / 300., 0.75);
    }
    inline double chiHeIICII(double T)
    {
        return 5e-20 * T * T * exp(-7e-6 * T) * exp(-73003. / T);
    }  // 6.29 eV

    // -- N --
    inline double xiHeINIII(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 1e7) ? 1000. : 1e-4 * T);
        if (T4 < 4.) return 4.84e-10 * pow(T4, 0.92) * (1. + 2.37 * exp(-10.2 * T4));
        return 3.17e-9 * pow(T4, 0.2) * (1. - 0.72 * exp(-0.0481 * T4));
    }
    inline double xiHeINIV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 1e7) ? 1000. : 1e-4 * T);
        return 2.05e-9 * pow(T4, 0.23) * (1. - 0.72 * exp(-0.19 * T4));
    }
    inline double xiHeINV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 1e7) ? 1000. : 1e-4 * T);
        if (T4 < 9.) return 1.26e-11 * pow(T4, 1.55) * (1. + 11.2 * exp(-7.82 * T4));
        return 3.75e-10 * pow(T4, 0.54) * (1. - 0.82 * exp(-0.0207 * T4));
    }
    inline double chiHeIINII(double T)
    {
        return 3.7e-20 * T * T * exp(-6.3e-6 * T) * exp(-16710. / T);
    }  // 1.44 eV

    // -- O --
    inline double xiHeIOIII(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 1e7) ? 1000. : 1e-4 * T);
        if (T4 < 5.) return 7.1e-12 * pow(T4, 2.6) * (1. + 8.99 * exp(-0.78 * T4));
        return 6.21e-10 * pow(T4, 0.53) * (1. - 0.66 * exp(-0.0222 * T4));
    }
    inline double xiHeIOIV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 1e7) ? 1000. : 1e-4 * T);
        return 1.12e-9 * pow(T4, 0.42) * (1. - 0.71 * exp(-0.0198 * T4));
    }
    inline double xiHeIOV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 1e7) ? 1000. : 1e-4 * T);
        return 9.97e-10 * pow(T4, 0.4) * (1. - 0.46 * exp(-0.35 * T4));
    }
    inline double chiHeIIOI(double T)
    {
        double T4 = 1e-4 * T;
        double T5 = (T < 1e7) ? 1e-5 * T : 100.;
        return 4.991e-15 * pow(T4, 0.3794) * exp(-T4 / 112.1) + 2.780e-15 * pow(T4, -0.2163) * exp(T5 / 8.158);
    }

    // -- Ne --
    inline double xiHeINeIII(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 1e7) ? 1000. : 1e-4 * T);
        if (T4 < 0.5) return 1e-14;
        if (T4 < 0.8) return 8.48e-12 * pow(T4, 3.35) * (1. - 1.92 * exp(-1.5 * T4));
        return 2.52e-11 * pow(T4, 0.14) * (1. - 1.99 * exp(-0.91 * T4));
    }
    inline double xiHeINeIV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 1e7) ? 1000. : 1e-4 * T);
        if (T4 < 2.5) return 1e-14;
        if (T4 < 9.5) return 1.34e-13 * pow(T4, 2.33) * (1. - 2.55 * exp(-0.37 * T4));
        return 1e-10 * pow(T4, 0.24) * (1. - 1.09 * exp(-0.0247 * T4));
    }
    inline double xiHeINeV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 1e7) ? 1000. : 1e-4 * T);
        if (T4 < 5.) return 1.77e-09 * pow(T4, 0.14) * (1. + 0.0488 * exp(-3.35 * T4));
        return 2.67e-10 * pow(T4, 0.54) * (1. + 0.91 * exp(-0.0188 * T4));
    }

    // -- Mg --
    inline double xiHeIMgIV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 1e7) ? 1000. : 1e-4 * T);
        if (T4 < 0.7) return 1.88e-11 * pow(T4, 0.13) * (1. + 0.83 * exp(-4.94 * T4));
        return 3.41e-11 * pow(T4, 0.15) * (1. - 0.45 * exp(-0.0483 * T4));
    }
    inline double xiHeIMgV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 1e7) ? 1000. : 1e-4 * T);
        return 1.37e-09 * pow(T4, 0.21) * (1. - 0.59 * exp(-0.0594 * T4));
    }

    // -- Si --
    inline double xiHeISiIV(double T)
    {
        double T4 = (T < 100.) ? 0.01 : ((T > 1e6) ? 100. : 1e-4 * T);
        return 1.03e-9 * pow(T4, 0.6) * (1. - 0.61 * exp(-1.42 * T4));
    }
    inline double xiHeISiV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 5e5) ? 50. : 1e-4 * T);
        return 5.75e-10 * pow(T4, 0.93) * (1. + 1.33 * exp(-0.29 * T4));
    }
    inline double chiHeIISiI(double /*T*/)
    {
        return 1.3e-9;
    }
    inline double chiHeIISiII(double T)
    {
        return 1.5e-11 * pow(T, 0.25) * exp(-80182. / T);
    }  // 6.91 eV
    inline double chiHeIISiIII(double T)
    {
        return 1.15e-11 * sqrt(T) * exp(-103052. / T);
    }  // 8.88 eV

    // -- S --
    inline double xiHeISIV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 3.1e4) ? 3.1 : 1e-4 * T);
        return 3.58e-9 * pow(T4, 0.00777) * (1. - 0.94 * exp(-0.3 * T4));
    }
    inline double xiHeISV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 3.1e4) ? 3.1 : 1e-4 * T);
        return 7.44e-13 * pow(T4, 0.34) * (1. + 3.74 * exp(-5.18 * T4));
    }
    inline double chiHeIISII(double T)
    {
        return 4.4e-16 * pow(T, 1.2) * exp(-3.6e-6 * T) * exp(-106764. / T);
    }  // 9.2 eV
    inline double chiHeIISIII(double T)
    {
        return 5.5e-18 * pow(T, 1.6) * exp(-4.6e-6 * T) * exp(-121853. / T);
    }  // 10.5 eV

    // -- Fe --
    inline double xiHeIFeIV(double /*T*/)
    {
        return 3.05e-11;
    }
    inline double xiHeIFeV(double T)
    {
        double T4 = (T < 1e3) ? 0.1 : ((T > 1e7) ? 1000. : 1e-4 * T);
        return 1.38e-9 * pow(T4, 0.43) * (1. + 1.12 * exp(-0.16 * T4));
    }

    /** Helium charge exchange recombination rate [cm^3/s] for a given ion (Verner index).
        Reaction: ion + HeI -> ion_lower + HeII.
        The ion index refers to the HIGHER stage (the one being recombined). */
    inline double xiHeI(int ion, double T)
    {
        switch (ion)
        {
            // C: CIV-CV (ions 6-7)
            case 6: return xiHeICIV(T);
            case 7: return xiHeICV(T);
            // N: NIII-NV (ions 11-13)
            case 11: return xiHeINIII(T);
            case 12: return xiHeINIV(T);
            case 13: return xiHeINV(T);
            // O: OIII-OV (ions 18-20)
            case 18: return xiHeIOIII(T);
            case 19: return xiHeIOIV(T);
            case 20: return xiHeIOV(T);
            // Ne: NeIII-NeV (ions 26-28)
            case 26: return xiHeINeIII(T);
            case 27: return xiHeINeIV(T);
            case 28: return xiHeINeV(T);
            // Mg: MgIV-MgV (ions 35-36)
            case 35: return xiHeIMgIV(T);
            case 36: return xiHeIMgV(T);
            // Si: SiIV-SiV (ions 45-46)
            case 45: return xiHeISiIV(T);
            case 46: return xiHeISiV(T);
            // S: SIV-SV (ions 56-57)
            case 56: return xiHeISIV(T);
            case 57: return xiHeISV(T);
            // Fe: FeIV-FeV (ions 61-62)
            case 61: return xiHeIFeIV(T);
            case 62: return xiHeIFeV(T);
            default: return 0.;
        }
    }

    /** Helium charge exchange ionization rate [cm^3/s] for a given ion (Verner index).
        Reaction: ion + HeII -> ion_higher + HeI.
        The ion index refers to the LOWER stage (the one being ionized). */
    inline double chiHeII(int ion, double T)
    {
        switch (ion)
        {
            // C: CI-CII (ions 3-4)
            case 3: return chiHeIICI(T);
            case 4: return chiHeIICII(T);
            // N: NII (ion 10)
            case 10: return chiHeIINII(T);
            // O: OI (ion 16)
            case 16: return chiHeIIOI(T);
            // Si: SiI-SiIII (ions 42-44)
            case 42: return chiHeIISiI(T);
            case 43: return chiHeIISiII(T);
            case 44: return chiHeIISiIII(T);
            // S: SII-SIII (ions 55-56)
            case 55: return chiHeIISII(T);
            case 56: return chiHeIISIII(T);
            default: return 0.;
        }
    }
}

//////////////////////////////////////////////////////////////////////

double PhotoIonizationRates::alphaTotal(int ion, double T)
{
    return ::alphaTotal(ion, T);
}
double PhotoIonizationRates::beta(int ion, double T)
{
    return ::beta(ion, T);
}
double PhotoIonizationRates::xiHI(int ion, double T)
{
    return ::xiHI(ion, T);
}
double PhotoIonizationRates::chiHII(int ion, double T)
{
    return ::chiHII(ion, T);
}
double PhotoIonizationRates::xiHeI(int ion, double T)
{
    return ::xiHeI(ion, T);
}
double PhotoIonizationRates::chiHeII(int ion, double T)
{
    return ::chiHeII(ion, T);
}

double PhotoIonizationRates::alphaBHII(double T)
{
    return ::alphaBHII(T);
}
double PhotoIonizationRates::alphaBHeII(double T)
{
    return ::alphaBHeII(T);
}
double PhotoIonizationRates::alphaBHeIII(double T)
{
    return ::alphaBHeIII(T);
}
double PhotoIonizationRates::alphaDRHeII(double T)
{
    return ::alphaDRHeII(T);
}

double PhotoIonizationRates::betaHI(double T)
{
    return ::betaHI(T);
}
double PhotoIonizationRates::betaHeI(double T)
{
    return ::betaHeI(T);
}
double PhotoIonizationRates::betaHeII(double T)
{
    return ::betaHeII(T);
}

double PhotoIonizationRates::xiHIHeII(double T)
{
    return ::xiHIHeII(T);
}
double PhotoIonizationRates::xiHIHeIII(double T)
{
    return ::xiHIHeIII(T);
}

//////////////////////////////////////////////////////////////////////
