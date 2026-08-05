/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FlatUniverseCosmology.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

double FlatUniverseCosmology::modelRedshift() const
{
    return redshift();
}

////////////////////////////////////////////////////////////////////

namespace
{
    // returns the comoving distance as defined in the class header, given the following parameters:
    //  z: redshift > 0
    //  h: reduced Hubble constant
    //  Om: matter density fraction
    double comovingDistance(double z, double h, double Om)
    {
        // front factor: c/H0
        double front = 10. * Constants::pc() * Constants::c() / h;

        // integral, using a fixed number of steps per redshift interval
        int n = max(2000, static_cast<int>(z * 10000));
        double dz = z / n;
        double sum = 0.;
        for (int i = 0; i != n; ++i)
        {
            double zp1 = (i + 0.5) * dz + 1.;
            double integrand = 1. / sqrt(Om * zp1 * zp1 * zp1 + (1. - Om));
            sum += integrand;
        }
        return front * z * sum / n;
    }
}

////////////////////////////////////////////////////////////////////

double FlatUniverseCosmology::angularDiameterDistance() const
{
    return comovingDistance(redshift(), reducedHubbleConstant(), matterDensityFraction()) / (1. + redshift());
}

////////////////////////////////////////////////////////////////////

double FlatUniverseCosmology::luminosityDistance() const
{
    return comovingDistance(redshift(), reducedHubbleConstant(), matterDensityFraction()) * (1. + redshift());
}

////////////////////////////////////////////////////////////////////

double FlatUniverseCosmology::relativeExpansionRate() const
{
    double zp1 = redshift() + 1;
    double H0 = reducedHubbleConstant() / (10. * Constants::pc());
    double Om = matterDensityFraction();
    return H0 * sqrt(Om * zp1 * zp1 * zp1 + (1. - Om));
}

////////////////////////////////////////////////////////////////////
