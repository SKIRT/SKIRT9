/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "VernerCrossSections.hpp"

//////////////////////////////////////////////////////////////////////

namespace
{
    using std::pow;
    using std::sqrt;

    // ---- hydrogen ----

    /** HI photoionization cross-section [cm^2] over [13.6, 50000] eV. */
    inline double sigmaHI(double E)
    {
        if (E < VernerCrossSections::HI_eV || E > 50000.) return 0.;
        double x = E / 0.4298;
        double xm1 = x - 1.;
        return 5.475e-14 * xm1 * xm1 * pow(x, -4.0185) / pow(1. + sqrt(x / 32.88), 2.963);
    }

    // ---- helium ----

    /** HeI photoionization cross-section [cm^2] over [24.59, 50000] eV. */
    inline double sigmaHeI(double E)
    {
        if (E < VernerCrossSections::HeI_eV || E > 50000.) return 0.;
        double x = E / 13.61 - 0.4434;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 4.562496);
        return 9.492e-16 * (xm1 * xm1 + 4.157521) * pow(y, -3.906) / pow(1. + sqrt(y / 1.469), 3.188);
    }

    /** HeII photoionization cross-section [cm^2] over [54.42, 50000] eV. */
    inline double sigmaHeII(double E)
    {
        if (E < VernerCrossSections::HeII_eV || E > 50000.) return 0.;
        double x = E / 1.72;
        double xm1 = x - 1.;
        return 1.369e-14 * xm1 * xm1 * pow(x, -4.0185) / pow(1. + sqrt(x / 32.88), 2.963);
    }

    // ---- carbon ----

    /** CI photoionization cross-section [cm^2] over [11.26, 291] eV. */
    inline double sigmaCI(double E)
    {
        if (E < VernerCrossSections::CI_eV || E > 291.) return 0.;
        double x = E / 2.144 - 1.133;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 2.582449);
        return 5.027e-16 * (xm1 * xm1 + 0.0083850649) * pow(y, -2.9495) / pow(1. + sqrt(y / 62.16), 5.101);
    }

    /** CII photoionization cross-section [cm^2] over [24.38, 307.6] eV. */
    inline double sigmaCII(double E)
    {
        if (E < VernerCrossSections::CII_eV || E > 307.6) return 0.;
        double x = E / 0.4058 - 49.29;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 10.458756);
        return 8.709e-18 * (xm1 * xm1 + 4.380649) * pow(y, -1.211) / pow(1. + sqrt(y / 126.1), 8.578);
    }

    /** CIII photoionization cross-section [cm^2] over [47.89, 328.9] eV. */
    inline double sigmaCIII(double E)
    {
        if (E < VernerCrossSections::CIII_eV || E > 328.9) return 0.;
        double x = E / 4.614 - 0.004378;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.0006390784);
        return 1.539e-14 * (xm1 * xm1 + 35.070084) * pow(y, 2.465) / pow(1. + sqrt(y / 1.737), 15.93);
    }

    /** CIV photoionization cross-section [cm^2] over [64.49, 352.2] eV. */
    inline double sigmaCIV(double E)
    {
        if (E < VernerCrossSections::CIV_eV || E > 352.2) return 0.;
        double x = E / 3.506;
        double xm1 = x - 1.;
        return 1.068e-16 * xm1 * xm1 * pow(x, -1.7715) / pow(1. + sqrt(x / 14.36), 7.457);
    }

    /** CV photoionization cross-section [cm^2] over [392.1, 50000] eV. */
    inline double sigmaCV(double E)
    {
        if (E < VernerCrossSections::CV_eV || E > 50000.) return 0.;
        double x = E / 46.24;
        double xm1 = x - 1.;
        return 2.344e-16 * xm1 * xm1 * pow(x, -4.2095) / pow(1. + sqrt(x / 21.83), 2.581);
    }

    /** CVI photoionization cross-section [cm^2] over [490, 50000] eV. */
    inline double sigmaCVI(double E)
    {
        if (E < VernerCrossSections::CVI_eV || E > 50000.) return 0.;
        double x = E / 15.48;
        double xm1 = x - 1.;
        return 1.521e-15 * xm1 * xm1 * pow(x, -4.0185) / pow(1. + sqrt(x / 32.88), 2.963);
    }

    // ---- nitrogen ----

    /** NI photoionization cross-section [cm^2] over [14.53, 404.8] eV. */
    inline double sigmaNI(double E)
    {
        if (E < VernerCrossSections::NI_eV || E > 404.8) return 0.;
        double x = E / 4.034 - 0.8598;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 5.405625);
        return 8.235e-16 * (xm1 * xm1 + 0.0082755409) * pow(y, -3.536) / pow(1. + sqrt(y / 80.33), 3.928);
    }

    /** NII photoionization cross-section [cm^2] over [29.6, 423.6] eV. */
    inline double sigmaNII(double E)
    {
        if (E < VernerCrossSections::NII_eV || E > 423.6) return 0.;
        double x = E / 0.06128 - 428.;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 412.09);
        return 1.944e-18 * (xm1 * xm1 + 108.7849) * pow(y, -1.1135) / pow(1. + sqrt(y / 816.3), 8.773);
    }

    /** NIII photoionization cross-section [cm^2] over [47.45, 447.3] eV. */
    inline double sigmaNIII(double E)
    {
        if (E < VernerCrossSections::NIII_eV || E > 447.3) return 0.;
        double x = E / 0.242 - 187.7;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 15.992001);
        return 9.375e-19 * (xm1 * xm1 + 3.4225) * pow(y, -0.922) / pow(1. + sqrt(y / 278.8), 9.156);
    }

    /** NIV photoionization cross-section [cm^2] over [77.47, 475.3] eV. */
    inline double sigmaNIV(double E)
    {
        if (E < VernerCrossSections::NIV_eV || E > 475.3) return 0.;
        double x = E / 5.494 - 0.006415;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.0003751969);
        return 1.69e-14 * (xm1 * xm1 + 62.473216) * pow(y, 3.03) / pow(1. + sqrt(y / 1.714), 17.06);
    }

    /** NV photoionization cross-section [cm^2] over [97.89, 504.3] eV. */
    inline double sigmaNV(double E)
    {
        if (E < VernerCrossSections::NV_eV || E > 504.3) return 0.;
        double x = E / 4.471;
        double xm1 = x - 1.;
        return 8.376e-17 * xm1 * xm1 * pow(x, -2.4985) / pow(1. + sqrt(x / 32.97), 6.003);
    }

    /** NVI photoionization cross-section [cm^2] over [552.1, 50000] eV. */
    inline double sigmaNVI(double E)
    {
        if (E < VernerCrossSections::NVI_eV || E > 50000.) return 0.;
        double x = E / 69.43;
        double xm1 = x - 1.;
        return 1.519e-16 * xm1 * xm1 * pow(x, -4.3425) / pow(1. + sqrt(x / 26.27), 2.315);
    }

    /** NVII photoionization cross-section [cm^2] over [667.1, 50000] eV. */
    inline double sigmaNVII(double E)
    {
        if (E < VernerCrossSections::NVII_eV || E > 50000.) return 0.;
        double x = E / 21.08;
        double xm1 = x - 1.;
        return 1.117e-15 * xm1 * xm1 * pow(x, -4.0185) / pow(1. + sqrt(x / 32.88), 2.963);
    }

    // ---- oxygen ----

    /** OI photoionization cross-section [cm^2] over [13.62, 538] eV. */
    inline double sigmaOI(double E)
    {
        if (E < VernerCrossSections::OI_eV || E > 538.) return 0.;
        double x = E / 1.24 - 8.698;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.01615441);
        return 1.745e-15 * (xm1 * xm1 + 0.0057592921) * pow(y, 3.32) / pow(1. + sqrt(y / 3.784), 17.64);
    }

    /** OII photoionization cross-section [cm^2] over [35.12, 558.1] eV. */
    inline double sigmaOII(double E)
    {
        if (E < VernerCrossSections::OII_eV || E > 558.1) return 0.;
        double x = E / 1.386 - 21.31;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.0002259009);
        return 5.967e-17 * (xm1 * xm1 + 0.0003740356) * pow(y, -1.0285) / pow(1. + sqrt(y / 31.75), 8.943);
    }

    /** OIII photoionization cross-section [cm^2] over [54.94, 584] eV. */
    inline double sigmaOIII(double E)
    {
        if (E < VernerCrossSections::OIII_eV || E > 584.) return 0.;
        double x = E / 0.1723 - 0.003839;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.20875761);
        return 6.753e-16 * (xm1 * xm1 + 0.01418481) * pow(y, -2.089) / pow(1. + sqrt(y / 385.2), 6.822);
    }

    /** OIV photoionization cross-section [cm^2] over [77.41, 614.4] eV. */
    inline double sigmaOIV(double E)
    {
        if (E < VernerCrossSections::OIV_eV || E > 614.4) return 0.;
        double x = E / 0.2044 - 332.8;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 1836.1225);
        return 8.659e-19 * (xm1 * xm1 + 9.878449) * pow(y, -1.1075) / pow(1. + sqrt(y / 493.1), 8.785);
    }

    /** OV photoionization cross-section [cm^2] over [113.9, 649.1] eV. */
    inline double sigmaOV(double E)
    {
        if (E < VernerCrossSections::OV_eV || E > 649.1) return 0.;
        double x = E / 2.854 - 0.03036;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.0030846916);
        return 1.642e-14 * (xm1 * xm1 + 804.2896) * pow(y, 7.735) / pow(1. + sqrt(y / 1.792), 26.47);
    }

    /** OVI photoionization cross-section [cm^2] over [138.1, 683.7] eV. */
    inline double sigmaOVI(double E)
    {
        if (E < VernerCrossSections::OVI_eV || E > 683.7) return 0.;
        double x = E / 7.824;
        double xm1 = x - 1.;
        return 6.864e-17 * xm1 * xm1 * pow(x, -2.7525) / pow(1. + sqrt(x / 32.1), 5.495);
    }

    /** OVII photoionization cross-section [cm^2] over [739.3, 50000] eV. */
    inline double sigmaOVII(double E)
    {
        if (E < VernerCrossSections::OVII_eV || E > 50000.) return 0.;
        double x = E / 87.09;
        double xm1 = x - 1.;
        return 1.329e-16 * xm1 * xm1 * pow(x, -4.332) / pow(1. + sqrt(x / 25.35), 2.336);
    }

    /** OVIII photoionization cross-section [cm^2] over [871.4, 50000] eV. */
    inline double sigmaOVIII(double E)
    {
        if (E < VernerCrossSections::OVIII_eV || E > 50000.) return 0.;
        double x = E / 27.54;
        double xm1 = x - 1.;
        return 8.554e-16 * xm1 * xm1 * pow(x, -4.0185) / pow(1. + sqrt(x / 32.88), 2.963);
    }

    // ---- neon ----

    /** NeI photoionization cross-section [cm^2] over [21.56, 870.1] eV. */
    inline double sigmaNeI(double E)
    {
        if (E < VernerCrossSections::NeI_eV || E > 870.1) return 0.;
        double x = E / 4.87 - 0.04236;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 34.492129);
        return 4.287e-15 * (xm1 * xm1 + 0.05924356) * pow(y, -1.3225) / pow(1. + sqrt(y / 5.798), 8.355);
    }

    /** NeII photoionization cross-section [cm^2] over [40.96, 883.1] eV. */
    inline double sigmaNeII(double E)
    {
        if (E < VernerCrossSections::NeII_eV || E > 883.1) return 0.;
        double x = E / 12.47 - 1.52;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.01175056);
        return 1.583e-15 * (xm1 * xm1 + 0.0043007364) * pow(y, -1.595) / pow(1. + sqrt(y / 3.935), 7.81);
    }

    /** NeIII photoionization cross-section [cm^2] over [63.46, 913.1] eV. */
    inline double sigmaNeIII(double E)
    {
        if (E < VernerCrossSections::NeIII_eV || E > 913.1) return 0.;
        double x = E / 0.7753 - 76.54;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 4.092529);
        return 5.708e-18 * (xm1 * xm1 + 0.21464689) * pow(y, -0.475) / pow(1. + sqrt(y / 67.25), 10.05);
    }

    /** NeIV photoionization cross-section [cm^2] over [97.12, 948] eV. */
    inline double sigmaNeIV(double E)
    {
        if (E < VernerCrossSections::NeIV_eV || E > 948.) return 0.;
        double x = E / 5.566 - 5.149;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 44.715969);
        return 1.685e-15 * (xm1 * xm1 + 6.87241e-05) * pow(y, -3.972) / pow(1. + sqrt(y / 640.9), 3.056);
    }

    /** NeV photoionization cross-section [cm^2] over [126.2, 987.3] eV. */
    inline double sigmaNeV(double E)
    {
        if (E < VernerCrossSections::NeV_eV || E > 987.3) return 0.;
        double x = E / 1.248 - 91.69;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.13704804);
        return 2.43e-18 * (xm1 * xm1 + 0.46991025) * pow(y, -1.0005) / pow(1. + sqrt(y / 106.6), 8.999);
    }

    /** NeVI photoionization cross-section [cm^2] over [157.9, 1031] eV. */
    inline double sigmaNeVI(double E)
    {
        if (E < VernerCrossSections::NeVI_eV || E > 1031.) return 0.;
        double x = E / 1.499 - 104.2;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 2.059225);
        return 9.854e-19 * (xm1 * xm1 + 2.742336) * pow(y, -1.082) / pow(1. + sqrt(y / 135.), 8.836);
    }

    /** NeVII photoionization cross-section [cm^2] over [207.3, 1078] eV. */
    inline double sigmaNeVII(double E)
    {
        if (E < VernerCrossSections::NeVII_eV || E > 1078.) return 0.;
        double x = E / 4.888 - 0.02536;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.0019509889);
        return 1.198e-14 * (xm1 * xm1 + 790.1721) * pow(y, 7.25) / pow(1. + sqrt(y / 1.788), 25.5);
    }

    /** NeVIII photoionization cross-section [cm^2] over [239.1, 1125] eV. */
    inline double sigmaNeVIII(double E)
    {
        if (E < VernerCrossSections::NeVIII_eV || E > 1125.) return 0.;
        double x = E / 10.03;
        double xm1 = x - 1.;
        return 5.631e-17 * xm1 * xm1 * pow(x, -2.7075) / pow(1. + sqrt(x / 36.28), 5.585);
    }

    // ---- magnesium (10 ions, MgI-MgX) ----

    /** MgI photoionization cross-section [cm^2] over [7.646, 54.9] eV. */
    inline double sigmaMgI(double E)
    {
        if (E < VernerCrossSections::MgI_eV || E > 54.9) return 0.;
        double x = E / 11.97;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0);
        return 1.372e-10 * (xm1 * xm1 + 0.07868025) * pow(y, 2.37) / pow(1. + sqrt(y / 0.2228), 15.74);
    }

    /** MgII photoionization cross-section [cm^2] over [15.04, 65.69] eV. */
    inline double sigmaMgII(double E)
    {
        if (E < VernerCrossSections::MgII_eV || E > 65.69) return 0.;
        double x = E / 8.139;
        double xm1 = x - 1.;
        return 3.278e-18 * xm1 * xm1 * pow(x, -3.695) / pow(1. + sqrt(x / 4.341e7), 3.61);
    }

    /** MgIII photoionization cross-section [cm^2] over [80.14, 1317] eV. */
    inline double sigmaMgIII(double E)
    {
        if (E < VernerCrossSections::MgIII_eV || E > 1317.) return 0.;
        double x = E / 10.86 - 4.86;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 13.853284);
        return 5.377e-16 * (xm1 * xm1 + 6.780816) * pow(y, -1.9415) / pow(1. + sqrt(y / 9.779), 7.117);
    }

    /** MgIV photoionization cross-section [cm^2] over [109.3, 1356] eV. */
    inline double sigmaMgIV(double E)
    {
        if (E < VernerCrossSections::MgIV_eV || E > 1356.) return 0.;
        double x = E / 29.12 - 0.9402;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.01288225);
        return 1.394e-15 * (xm1 * xm1 + 0.0018714276) * pow(y, -2.2565) / pow(1. + sqrt(y / 2.895), 6.487);
    }

    /** MgV photoionization cross-section [cm^2] over [141.3, 1400] eV. */
    inline double sigmaMgV(double E)
    {
        if (E < VernerCrossSections::MgV_eV || E > 1400.) return 0.;
        double x = E / 0.9762 - 127.6;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 15.832441);
        return 1.728e-18 * (xm1 * xm1 + 0.654481) * pow(y, -0.47) / pow(1. + sqrt(y / 91.84), 10.06);
    }

    /** MgVI photoionization cross-section [cm^2] over [186.5, 1449] eV. */
    inline double sigmaMgVI(double E)
    {
        if (E < VernerCrossSections::MgVI_eV || E > 1449.) return 0.;
        double x = E / 1.711 - 100.7;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 2.989441);
        return 2.185e-18 * (xm1 * xm1 + 0.40005625) * pow(y, -0.899) / pow(1. + sqrt(y / 93.5), 9.202);
    }

    /** MgVII photoionization cross-section [cm^2] over [224.9, 1501] eV. */
    inline double sigmaMgVII(double E)
    {
        if (E < VernerCrossSections::MgVII_eV || E > 1501.) return 0.;
        double x = E / 3.57 - 54.52;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 4.318084);
        return 3.104e-18 * (xm1 * xm1 + 2.022084) * pow(y, -1.0715) / pow(1. + sqrt(y / 60.6), 8.857);
    }

    /** MgVIII photoionization cross-section [cm^2] over [266, 1558] eV. */
    inline double sigmaMgVIII(double E)
    {
        if (E < VernerCrossSections::MgVIII_eV || E > 1558.) return 0.;
        double x = E / 0.4884 - 534.8;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 1.5976009e-05);
        return 6.344e-20 * (xm1 * xm1 + 0.44435556) * pow(y, -0.8075) / pow(1. + sqrt(y / 508.5), 9.385);
    }

    /** MgIX photoionization cross-section [cm^2] over [328.2, 1618] eV. */
    inline double sigmaMgIX(double E)
    {
        if (E < VernerCrossSections::MgIX_eV || E > 1618.) return 0.;
        double x = E / 34.82 - 5.444;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.0062694724);
        return 9.008e-16 * (xm1 * xm1 + 7.568001) * pow(y, 1.72) / pow(1. + sqrt(y / 1.823), 14.44);
    }

    /** MgX photoionization cross-section [cm^2] over [367.5, 1675] eV. */
    inline double sigmaMgX(double E)
    {
        if (E < VernerCrossSections::MgX_eV || E > 1675.) return 0.;
        double x = E / 14.52;
        double xm1 = x - 1.;
        return 4.427e-17 * xm1 * xm1 * pow(x, -2.77) / pow(1. + sqrt(x / 38.26), 5.46);
    }

    // ---- silicon (12 ions, SiI-SiXII) ----

    /** SiI photoionization cross-section [cm^2] over [8.152, 106] eV. */
    inline double sigmaSiI(double E)
    {
        if (E < VernerCrossSections::SiI_eV || E > 106.) return 0.;
        double x = E / 23.17 - 1.672e-05;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.17698849);
        return 2.506e-17 * (xm1 * xm1 + 0.08048569) * pow(y, -3.727) / pow(1. + sqrt(y / 20.57), 3.546);
    }

    /** SiII photoionization cross-section [cm^2] over [16.35, 118.6] eV. */
    inline double sigmaSiII(double E)
    {
        if (E < VernerCrossSections::SiII_eV || E > 118.6) return 0.;
        double x = E / 2.556 - 6.634;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.01617984);
        return 4.14e-18 * (xm1 * xm1 + 2.4649) * pow(y, 0.455) / pow(1. + sqrt(y / 13.37), 11.91);
    }

    /** SiIII photoionization cross-section [cm^2] over [33.49, 131.1] eV. */
    inline double sigmaSiIII(double E)
    {
        if (E < VernerCrossSections::SiIII_eV || E > 131.1) return 0.;
        double x = E / 0.1659 - 96.13;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.41499364);
        return 5.79e-22 * (xm1 * xm1 + 0.74407876) * pow(y, 1.18) / pow(1. + sqrt(y / 147.4), 13.36);
    }

    /** SiIV photoionization cross-section [cm^2] over [45.14, 146.6] eV. */
    inline double sigmaSiIV(double E)
    {
        if (E < VernerCrossSections::SiIV_eV || E > 146.6) return 0.;
        double x = E / 12.88;
        double xm1 = x - 1.;
        return 6.083e-18 * xm1 * xm1 * pow(x, -3.8235) / pow(1. + sqrt(x / 1.356e6), 3.353);
    }

    /** SiV photoionization cross-section [cm^2] over [166.8, 1887] eV. */
    inline double sigmaSiV(double E)
    {
        if (E < VernerCrossSections::SiV_eV || E > 1887.) return 0.;
        double x = E / 0.7761 - 200.9;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 20.584369);
        return 8.863e-19 * (xm1 * xm1 + 1.697809) * pow(y, -0.51) / pow(1. + sqrt(y / 154.1), 9.98);
    }

    /** SiVI photoionization cross-section [cm^2] over [205.1, 1946] eV. */
    inline double sigmaSiVI(double E)
    {
        if (E < VernerCrossSections::SiVI_eV || E > 1946.) return 0.;
        double x = E / 63.05 - 1.115;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.0064818601);
        return 7.293e-17 * (xm1 * xm1 + 8.934121e-06) * pow(y, -4.3) / pow(1. + sqrt(y / 155.8), 2.4);
    }

    /** SiVII photoionization cross-section [cm^2] over [246.5, 2001] eV. */
    inline double sigmaSiVII(double E)
    {
        if (E < VernerCrossSections::SiVII_eV || E > 2001.) return 0.;
        double x = E / 0.3277 - 0.01149;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.40908816);
        return 6.68e-20 * (xm1 * xm1 + 10.7584) * pow(y, 2.53) / pow(1. + sqrt(y / 41.32), 16.06);
    }

    /** SiVIII photoionization cross-section [cm^2] over [303.2, 2058] eV. */
    inline double sigmaSiVIII(double E)
    {
        if (E < VernerCrossSections::SiVIII_eV || E > 2058.) return 0.;
        double x = E / 0.7655 - 385.;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.0080982001);
        return 3.477e-19 * (xm1 * xm1 + 2.178576e-06) * pow(y, -1.007) / pow(1. + sqrt(y / 373.3), 8.986);
    }

    /** SiIX photoionization cross-section [cm^2] over [351.1, 2125] eV. */
    inline double sigmaSiIX(double E)
    {
        if (E < VernerCrossSections::SiIX_eV || E > 2125.) return 0.;
        double x = E / 0.3343 - 1036.;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.08620096);
        return 1.465e-19 * (xm1 * xm1 + 2.709316) * pow(y, -1.2485) / pow(1. + sqrt(y / 1404.), 8.503);
    }

    /** SiX photoionization cross-section [cm^2] over [401.4, 2194] eV. */
    inline double sigmaSiX(double E)
    {
        if (E < VernerCrossSections::SiX_eV || E > 2194.) return 0.;
        double x = E / 0.8787 - 452.8;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 1.030225);
        return 1.95e-19 * (xm1 * xm1 + 0.20151121) * pow(y, -1.349) / pow(1. + sqrt(y / 746.1), 8.302);
    }

    /** SiXI photoionization cross-section [cm^2] over [476.1, 2268] eV. */
    inline double sigmaSiXI(double E)
    {
        if (E < VernerCrossSections::SiXI_eV || E > 2268.) return 0.;
        double x = E / 12.05 - 0.0199;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.0001014049);
        return 1.992e-14 * (xm1 * xm1 + 572.1664) * pow(y, 6.625) / pow(1. + sqrt(y / 1.582), 24.25);
    }

    /** SiXII photoionization cross-section [cm^2] over [523.5, 2336] eV. */
    inline double sigmaSiXII(double E)
    {
        if (E < VernerCrossSections::SiXII_eV || E > 2336.) return 0.;
        double x = E / 35.6;
        double xm1 = x - 1.;
        return 2.539e-17 * xm1 * xm1 * pow(x, -3.136) / pow(1. + sqrt(x / 33.07), 4.728);
    }

    // ---- sulfur (first 4 ions, matching CHIANTI coverage) ----

    /** SI photoionization cross-section [cm^2] over [10.36, 170] eV. */
    inline double sigmaSI(double E)
    {
        if (E < VernerCrossSections::SI_eV || E > 170.) return 0.;
        double x = E / 18.08 - 0.9935;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.06180196);
        return 4.564e-14 * (xm1 * xm1 + 0.40768225) * pow(y, 1.305) / pow(1. + sqrt(y), 13.61);
    }

    /** SII photoionization cross-section [cm^2] over [23.33, 184.6] eV. */
    inline double sigmaSII(double E)
    {
        if (E < VernerCrossSections::SII_eV || E > 184.6) return 0.;
        double x = E / 8.787 - 2.782;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.03196944);
        return 3.136e-16 * (xm1 * xm1 + 0.54081316) * pow(y, 0.905) / pow(1. + sqrt(y / 3.442), 12.81);
    }

    /** SIII photoionization cross-section [cm^2] over [34.83, 199.5] eV. */
    inline double sigmaSIII(double E)
    {
        if (E < VernerCrossSections::SIII_eV || E > 199.5) return 0.;
        double x = E / 2.027 - 15.68;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 88.755241);
        return 6.666e-18 * (xm1 * xm1 + 16.883881) * pow(y, -1.1945) / pow(1. + sqrt(y / 54.54), 8.611);
    }

    /** SIV photoionization cross-section [cm^2] over [47.31, 216.4] eV. */
    inline double sigmaSIV(double E)
    {
        if (E < VernerCrossSections::SIV_eV || E > 216.4) return 0.;
        double x = E / 2.173 - 19.75;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 11.296321);
        return 2.606e-18 * (xm1 * xm1 + 3.470769) * pow(y, -1.1725) / pow(1. + sqrt(y / 66.41), 8.655);
    }

    // ---- iron (16 ionization stages, matching CHIANTI coverage) ----

    /** FeI photoionization cross-section [cm^2] over [7.902, 66] eV. */
    inline double sigmaFeI(double E)
    {
        if (E < VernerCrossSections::FeI_eV || E > 66.) return 0.;
        double x = E / 0.05461 - 138.2;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.06155361);
        return 3.062e-19 * (xm1 * xm1 + 428.0761) * pow(y, -1.5385) / pow(1. + sqrt(y / 2.671e7), 7.923);
    }

    /** FeII photoionization cross-section [cm^2] over [16.19, 76.17] eV. */
    inline double sigmaFeII(double E)
    {
        if (E < VernerCrossSections::FeII_eV || E > 76.17) return 0.;
        double x = E / 0.1761 - 92.72;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 11556.25);
        return 4.365e-15 * (xm1 * xm1 + 130.1881) * pow(y, -2.898) / pow(1. + sqrt(y / 6298.), 5.204);
    }

    /** FeIII photoionization cross-section [cm^2] over [30.65, 87.05] eV. */
    inline double sigmaFeIII(double E)
    {
        if (E < VernerCrossSections::FeIII_eV || E > 87.05) return 0.;
        double x = E / 0.1698 - 176.;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 341.1409);
        return 6.107e-18 * (xm1 * xm1 + 75.655204) * pow(y, -1.4725) / pow(1. + sqrt(y / 1555.), 8.055);
    }

    /** FeIV photoionization cross-section [cm^2] over [54.8, 106.7] eV. */
    inline double sigmaFeIV(double E)
    {
        if (E < VernerCrossSections::FeIV_eV || E > 106.7) return 0.;
        double x = E / 25.44;
        double xm1 = x - 1.;
        return 3.653e-16 * (xm1 * xm1 + 0.31382404) * pow(x, -2.231) / pow(1. + sqrt(x / 8.913), 6.538);
    }

    /** FeV photoionization cross-section [cm^2] over [75.01, 128.8] eV. */
    inline double sigmaFeV(double E)
    {
        if (E < VernerCrossSections::FeV_eV || E > 128.8) return 0.;
        double x = E / 0.7256 - 88.71;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.00278784);
        return 1.523e-21 * (xm1 * xm1 + 2564.4096) * pow(y, 3.335) / pow(1. + sqrt(y / 37.36), 17.67);
    }

    /** FeVI photoionization cross-section [cm^2] over [99.06, 152.7] eV. */
    inline double sigmaFeVI(double E)
    {
        if (E < VernerCrossSections::FeVI_eV || E > 152.7) return 0.;
        double x = E / 2.656 - 33.61;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 1.4010049e-05);
        return 5.259e-19 * (xm1 * xm1 + 242.7364) * pow(y, 2.66) / pow(1. + sqrt(y / 14.5), 16.32);
    }

    /** FeVII photoionization cross-section [cm^2] over [125, 178.3] eV. */
    inline double sigmaFeVII(double E)
    {
        if (E < VernerCrossSections::FeVII_eV || E > 178.3) return 0.;
        double x = E / 5.059 - 0.4546;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 719.8489);
        return 2.42e-14 * (xm1 * xm1 + 6.330256e-06) * pow(y, -4.313) / pow(1. + sqrt(y / 48500.), 2.374);
    }

    /** FeVIII photoionization cross-section [cm^2] over [151.1, 205.5] eV. */
    inline double sigmaFeVIII(double E)
    {
        if (E < VernerCrossSections::FeVIII_eV || E > 205.5) return 0.;
        double x = E / 0.07098 - 2542.;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 218275.84);
        return 1.979e-17 * (xm1 * xm1 + 46569.64) * pow(y, -2.125) / pow(1. + sqrt(y / 17450.), 6.75);
    }

    /** FeIX photoionization cross-section [cm^2] over [233.6, 921.1] eV. */
    inline double sigmaFeIX(double E)
    {
        if (E < VernerCrossSections::FeIX_eV || E > 921.1) return 0.;
        double x = E / 6.741 - 24.94;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 68.079001);
        return 2.687e-17 * (xm1 * xm1 + 5.697769e-08) * pow(y, -2.355) / pow(1. + sqrt(y / 180.7), 6.29);
    }

    /** FeX photoionization cross-section [cm^2] over [262.1, 959] eV. */
    inline double sigmaFeX(double E)
    {
        if (E < VernerCrossSections::FeX_eV || E > 959.) return 0.;
        double x = E / 68.86 - 1.19e-05;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 4.31649e-05);
        return 6.47e-17 * (xm1 * xm1 + 7.717284e-08) * pow(y, -3.4445) / pow(1. + sqrt(y / 20.62), 4.111);
    }

    /** FeXI photoionization cross-section [cm^2] over [290.2, 998.3] eV. */
    inline double sigmaFeXI(double E)
    {
        if (E < VernerCrossSections::FeXI_eV || E > 998.3) return 0.;
        double x = E / 8.284 - 29.71;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.272484);
        return 3.281e-18 * (xm1 * xm1 + 0.10751841) * pow(y, -1.2145) / pow(1. + sqrt(y / 53.6), 8.571);
    }

    /** FeXII photoionization cross-section [cm^2] over [330.8, 1039] eV. */
    inline double sigmaFeXII(double E)
    {
        if (E < VernerCrossSections::FeXII_eV || E > 1039.) return 0.;
        double x = E / 6.295 - 46.71;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 0.02030625);
        return 1.738e-18 * (xm1 * xm1 + 0.09585216) * pow(y, -1.4815) / pow(1. + sqrt(y / 113.), 8.037);
    }

    /** FeXIII photoionization cross-section [cm^2] over [361, 1081] eV. */
    inline double sigmaFeXIII(double E)
    {
        if (E < VernerCrossSections::FeXIII_eV || E > 1081.) return 0.;
        double x = E / 0.1317 - 2170.;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 4.6949904e-05);
        return 2.791e-21 * (xm1 * xm1 + 0.48135844) * pow(y, -0.6045) / pow(1. + sqrt(y / 2487.), 9.791);
    }

    /** FeXIV photoionization cross-section [cm^2] over [392.2, 1125] eV. */
    inline double sigmaFeXIV(double E)
    {
        if (E < VernerCrossSections::FeXIV_eV || E > 1125.) return 0.;
        double x = E / 0.8509 - 450.5;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 6.270016);
        return 1.454e-19 * (xm1 * xm1 + 0.24373969) * pow(y, -1.467) / pow(1. + sqrt(y / 1239.), 8.066);
    }

    /** FeXV photoionization cross-section [cm^2] over [457, 1181] eV. */
    inline double sigmaFeXV(double E)
    {
        if (E < VernerCrossSections::FeXV_eV || E > 1181.) return 0.;
        double x = E / 0.05555 - 0.0002706;
        double xm1 = x - 1.;
        double y = sqrt(x * x + 2.650384);
        return 2.108e-16 * (xm1 * xm1 + 3.553225e-06) * pow(y, -2.4835) / pow(1. + sqrt(y / 20450.), 6.033);
    }

    /** FeXVI photoionization cross-section [cm^2] over [489.3, 1216] eV. */
    inline double sigmaFeXVI(double E)
    {
        if (E < VernerCrossSections::FeXVI_eV || E > 1216.) return 0.;
        double x = E / 28.73;
        double xm1 = x - 1.;
        return 1.207e-17 * xm1 * xm1 * pow(x, -3.577) / pow(1. + sqrt(x / 515.), 3.846);
    }
}

//////////////////////////////////////////////////////////////////////

double VernerCrossSections::sigmaHI(double E)
{
    return ::sigmaHI(E);
}

//////////////////////////////////////////////////////////////////////

double VernerCrossSections::sigmaHeI(double E)
{
    return ::sigmaHeI(E);
}

//////////////////////////////////////////////////////////////////////

double VernerCrossSections::sigmaHeII(double E)
{
    return ::sigmaHeII(E);
}

//////////////////////////////////////////////////////////////////////

double VernerCrossSections::ionizationPotential(int ion)
{
    static constexpr double ip[] = {
        HI_eV,    HeI_eV,   HeII_eV, CI_eV,    CII_eV,    CIII_eV,   CIV_eV,    CV_eV,   CVI_eV,   NI_eV,     NII_eV,
        NIII_eV,  NIV_eV,   NV_eV,   NVI_eV,   NVII_eV,   OI_eV,     OII_eV,    OIII_eV, OIV_eV,   OV_eV,     OVI_eV,
        OVII_eV,  OVIII_eV, NeI_eV,  NeII_eV,  NeIII_eV,  NeIV_eV,   NeV_eV,    NeVI_eV, NeVII_eV, NeVIII_eV, MgI_eV,
        MgII_eV,  MgIII_eV, MgIV_eV, MgV_eV,   MgVI_eV,   MgVII_eV,  MgVIII_eV, MgIX_eV, MgX_eV,   SiI_eV,    SiII_eV,
        SiIII_eV, SiIV_eV,  SiV_eV,  SiVI_eV,  SiVII_eV,  SiVIII_eV, SiIX_eV,   SiX_eV,  SiXI_eV,  SiXII_eV,  SI_eV,
        SII_eV,   SIII_eV,  SIV_eV,  FeI_eV,   FeII_eV,   FeIII_eV,  FeIV_eV,   FeV_eV,  FeVI_eV,  FeVII_eV,  FeVIII_eV,
        FeIX_eV,  FeX_eV,   FeXI_eV, FeXII_eV, FeXIII_eV, FeXIV_eV,  FeXV_eV,   FeXVI_eV};
    return ip[ion];
}

//////////////////////////////////////////////////////////////////////

double VernerCrossSections::sigma(int ion, double E)
{
    switch (ion)
    {
        case 0: return ::sigmaHI(E);
        case 1: return ::sigmaHeI(E);
        case 2: return ::sigmaHeII(E);
        case 3: return ::sigmaCI(E);
        case 4: return ::sigmaCII(E);
        case 5: return ::sigmaCIII(E);
        case 6: return ::sigmaCIV(E);
        case 7: return ::sigmaCV(E);
        case 8: return ::sigmaCVI(E);
        case 9: return ::sigmaNI(E);
        case 10: return ::sigmaNII(E);
        case 11: return ::sigmaNIII(E);
        case 12: return ::sigmaNIV(E);
        case 13: return ::sigmaNV(E);
        case 14: return ::sigmaNVI(E);
        case 15: return ::sigmaNVII(E);
        case 16: return ::sigmaOI(E);
        case 17: return ::sigmaOII(E);
        case 18: return ::sigmaOIII(E);
        case 19: return ::sigmaOIV(E);
        case 20: return ::sigmaOV(E);
        case 21: return ::sigmaOVI(E);
        case 22: return ::sigmaOVII(E);
        case 23: return ::sigmaOVIII(E);
        case 24: return ::sigmaNeI(E);
        case 25: return ::sigmaNeII(E);
        case 26: return ::sigmaNeIII(E);
        case 27: return ::sigmaNeIV(E);
        case 28: return ::sigmaNeV(E);
        case 29: return ::sigmaNeVI(E);
        case 30: return ::sigmaNeVII(E);
        case 31: return ::sigmaNeVIII(E);
        case 32: return ::sigmaMgI(E);
        case 33: return ::sigmaMgII(E);
        case 34: return ::sigmaMgIII(E);
        case 35: return ::sigmaMgIV(E);
        case 36: return ::sigmaMgV(E);
        case 37: return ::sigmaMgVI(E);
        case 38: return ::sigmaMgVII(E);
        case 39: return ::sigmaMgVIII(E);
        case 40: return ::sigmaMgIX(E);
        case 41: return ::sigmaMgX(E);
        case 42: return ::sigmaSiI(E);
        case 43: return ::sigmaSiII(E);
        case 44: return ::sigmaSiIII(E);
        case 45: return ::sigmaSiIV(E);
        case 46: return ::sigmaSiV(E);
        case 47: return ::sigmaSiVI(E);
        case 48: return ::sigmaSiVII(E);
        case 49: return ::sigmaSiVIII(E);
        case 50: return ::sigmaSiIX(E);
        case 51: return ::sigmaSiX(E);
        case 52: return ::sigmaSiXI(E);
        case 53: return ::sigmaSiXII(E);
        case 54: return ::sigmaSI(E);
        case 55: return ::sigmaSII(E);
        case 56: return ::sigmaSIII(E);
        case 57: return ::sigmaSIV(E);
        case 58: return ::sigmaFeI(E);
        case 59: return ::sigmaFeII(E);
        case 60: return ::sigmaFeIII(E);
        case 61: return ::sigmaFeIV(E);
        case 62: return ::sigmaFeV(E);
        case 63: return ::sigmaFeVI(E);
        case 64: return ::sigmaFeVII(E);
        case 65: return ::sigmaFeVIII(E);
        case 66: return ::sigmaFeIX(E);
        case 67: return ::sigmaFeX(E);
        case 68: return ::sigmaFeXI(E);
        case 69: return ::sigmaFeXII(E);
        case 70: return ::sigmaFeXIII(E);
        case 71: return ::sigmaFeXIV(E);
        case 72: return ::sigmaFeXV(E);
        case 73: return ::sigmaFeXVI(E);
        default: return 0.;
    }
}

//////////////////////////////////////////////////////////////////////
