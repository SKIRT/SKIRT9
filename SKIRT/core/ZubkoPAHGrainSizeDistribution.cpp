/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ZubkoPAHGrainSizeDistribution.hpp"

//////////////////////////////////////////////////////////////////////

namespace
{
    // grain size range for PAH (in m)
    const double amin_pah = 0.00035e-6;
    const double amax_pah = 0.005e-6;

    // parameterized Zubko grain size distribution
    double dnda(double a, double A, double c0, double b0, double a1, double b1, double m1, double a2, double b2,
                double m2, double a3, double b3, double m3, double a4, double b4, double m4)
    {
        a *= 1e6;  // convert from m to micron
        double logg = c0 + b0 * log10(a) - b1 * pow(fabs(log10(a / a1)), m1) - b2 * pow(fabs(log10(a / a2)), m2)
                      - b3 * pow(fabs(a - a3), m3) - b4 * pow(fabs(a - a4), m4);
        return 1e6 * A * pow(10.0, logg);  // convert from 1/micron to 1/m
    }

    // grain size distribution for PAH
    double dnda_pah(double a)
    {
        const double A = 2.227433e-7;
        const double c0 = -8.02895;
        const double b0 = -3.45764;
        const double a1 = 1.0;  // in micron
        const double b1 = 1183.96;
        const double m1 = -8.20551;
        const double a2 = 1.0;  // not used
        const double b2 = 0.0;  // not used
        const double m2 = 0.0;  // not used
        const double a3 = -5.29496e-3;
        const double b3 = 1.0e24;
        const double m3 = 12.0146;
        const double a4 = 1.0;  // not used
        const double b4 = 0.0;  // not used
        const double m4 = 0.0;  // not used
        return dnda(a, A, c0, b0, a1, b1, m1, a2, b2, m2, a3, b3, m3, a4, b4, m4);
    }
}

////////////////////////////////////////////////////////////////////

ZubkoPAHGrainSizeDistribution::ZubkoPAHGrainSizeDistribution(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

double ZubkoPAHGrainSizeDistribution::amin() const
{
    return amin_pah;
}

////////////////////////////////////////////////////////////////////

double ZubkoPAHGrainSizeDistribution::amax() const
{
    return amax_pah;
}

////////////////////////////////////////////////////////////////////

double ZubkoPAHGrainSizeDistribution::dnda(double a) const
{
    return dnda_pah(a);
}

////////////////////////////////////////////////////////////////////
