/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ZubkoGraphiteGrainSizeDistribution.hpp"

//////////////////////////////////////////////////////////////////////

namespace
{
    // grain size range for graphite (in m)
    const double amin_gra = 0.00035e-6;
    const double amax_gra = 0.33e-6;

    // parameterized Zubko grain size distribution
    double dnda(double a, double A, double c0, double b0, double a1, double b1, double m1, double a2, double b2,
                double m2, double a3, double b3, double m3, double a4, double b4, double m4)
    {
        a *= 1e6;  // convert from m to micron
        double logg = c0 + b0 * log10(a) - b1 * pow(fabs(log10(a / a1)), m1) - b2 * pow(fabs(log10(a / a2)), m2)
                      - b3 * pow(fabs(a - a3), m3) - b4 * pow(fabs(a - a4), m4);
        return 1e6 * A * pow(10.0, logg);  // convert from 1/micron to 1/m
    }

    // grain size distribution for graphite
    double dnda_gra(double a)
    {
        const double A = 1.905816e-7;
        const double c0 = -9.86;
        const double b0 = -5.02082;
        const double a1 = 0.415861;  // in micron
        const double b1 = 5.81215e-3;
        const double m1 = 4.63229;
        const double a2 = 1.0;  // not used
        const double b2 = 0.0;  // not used
        const double m2 = 0.0;  // not used
        const double a3 = 0.160344;
        const double b3 = 1125.02;
        const double m3 = 3.69897;
        const double a4 = 0.160501;
        const double b4 = 1126.02;
        const double m4 = 3.69967;
        return dnda(a, A, c0, b0, a1, b1, m1, a2, b2, m2, a3, b3, m3, a4, b4, m4);
    }
}

////////////////////////////////////////////////////////////////////

ZubkoGraphiteGrainSizeDistribution::ZubkoGraphiteGrainSizeDistribution(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

double ZubkoGraphiteGrainSizeDistribution::amin() const
{
    return amin_gra;
}

////////////////////////////////////////////////////////////////////

double ZubkoGraphiteGrainSizeDistribution::amax() const
{
    return amax_gra;
}

////////////////////////////////////////////////////////////////////

double ZubkoGraphiteGrainSizeDistribution::dnda(double a) const
{
    return dnda_gra(a);
}

////////////////////////////////////////////////////////////////////
