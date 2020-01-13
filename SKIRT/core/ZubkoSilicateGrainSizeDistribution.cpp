/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ZubkoSilicateGrainSizeDistribution.hpp"

//////////////////////////////////////////////////////////////////////

namespace
{
    // grain size range for silicate (in m)
    const double amin_sil = 0.00035e-6;
    const double amax_sil = 0.37e-6;

    // parameterized Zubko grain size distribution
    double dnda(double a, double A, double c0, double b0, double a1, double b1, double m1, double a2, double b2,
                double m2, double a3, double b3, double m3, double a4, double b4, double m4)
    {
        a *= 1e6;  // convert from m to micron
        double logg = c0 + b0 * log10(a) - b1 * pow(fabs(log10(a / a1)), m1) - b2 * pow(fabs(log10(a / a2)), m2)
                      - b3 * pow(fabs(a - a3), m3) - b4 * pow(fabs(a - a4), m4);
        return 1e6 * A * pow(10.0, logg);  // convert from 1/micron to 1/m
    }

    // grain size distribution for silicate
    double dnda_sil(double a)
    {
        const double A = 1.471288e-7;
        const double c0 = -8.47091;
        const double b0 = -3.68708;
        const double a1 = 7.64943e-3;  // in micron
        const double b1 = 2.37316e-5;
        const double m1 = 22.5489;
        const double a2 = 1.0;  // not used
        const double b2 = 0.0;  // not used
        const double m2 = 0.0;  // not used
        const double a3 = 0.480229;
        const double b3 = 2961.28;
        const double m3 = 12.1717;
        const double a4 = 1.0;  // not used
        const double b4 = 0.0;  // not used
        const double m4 = 0.0;  // not used
        return dnda(a, A, c0, b0, a1, b1, m1, a2, b2, m2, a3, b3, m3, a4, b4, m4);
    }
}

////////////////////////////////////////////////////////////////////

ZubkoSilicateGrainSizeDistribution::ZubkoSilicateGrainSizeDistribution(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

double ZubkoSilicateGrainSizeDistribution::amin() const
{
    return amin_sil;
}

////////////////////////////////////////////////////////////////////

double ZubkoSilicateGrainSizeDistribution::amax() const
{
    return amax_sil;
}

////////////////////////////////////////////////////////////////////

double ZubkoSilicateGrainSizeDistribution::dnda(double a) const
{
    return dnda_sil(a);
}

////////////////////////////////////////////////////////////////////
