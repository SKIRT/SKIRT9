/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SingleGrainDustMix.hpp"
#include "StoredTable.hpp"

////////////////////////////////////////////////////////////////////

double SingleGrainDustMix::getOpticalProperties(const Array& lambdav, const Array& thetav, Array& sigmaabsv,
                                                Array& sigmascav, Array& asymmparv, Table<2>& S11vv, Table<2>& S12vv,
                                                Table<2>& S33vv, Table<2>& S34vv, ArrayTable<2>& /*sigmaabsvv*/,
                                                ArrayTable<2>& /*sigmaabspolvv*/)
{
    // open the stored table file containing the basic optical properties
    string opticalPropsName = resourceNameForOpticalProps();
    StoredTable<1> sigmaabs(this, opticalPropsName, "lambda(m)", "sigmaabs(m2/H)");
    StoredTable<1> sigmasca(this, opticalPropsName, "lambda(m)", "sigmasca(m2/H)");
    StoredTable<1> asymmpar(this, opticalPropsName, "lambda(m)", "g(1)");

    // log warning if the simulation wavelength range extends beyond the optical property range
    informAvailableWavelengthRange(sigmaabs.axisRange<0>());

    // retrieve the optical properties on the requested wavelength grid
    int numLambda = lambdav.size();
    for (int ell = 0; ell != numLambda; ++ell)
    {
        double lambda = lambdav[ell];
        sigmaabsv[ell] = sigmaabs(lambda);
        sigmascav[ell] = sigmasca(lambda);
        asymmparv[ell] = asymmpar(lambda);
    }

    // get mu value for some arbitrary wavelength
    double mu = StoredTable<1>(this, opticalPropsName, "lambda(m)", "mu(kg/H)")[1.];

    // get the scattering mode advertised by this dust mix
    auto mode = scatteringMode();
    if (mode == ScatteringMode::MaterialPhaseFunction || mode == ScatteringMode::SphericalPolarization
        || mode == ScatteringMode::SpheroidalPolarization)
    {
        // open the stored table file containing the Mueller matrix coefficients
        string muellerName = resourceNameForMuellerMatrix();
        StoredTable<2> S11, S12, S33, S34;
        S11.open(this, muellerName, "lambda(m),theta(rad)", "S11(1)");
        if (mode == ScatteringMode::SphericalPolarization || mode == ScatteringMode::SpheroidalPolarization)
        {
            S12.open(this, muellerName, "lambda(m),theta(rad)", "S12(1)");
            S33.open(this, muellerName, "lambda(m),theta(rad)", "S33(1)");
            S34.open(this, muellerName, "lambda(m),theta(rad)", "S34(1)");
        }

        // retrieve the Mueller matrix coefficients on the requested wavelength and scattering angle grid
        int numTheta = thetav.size();
        for (int ell = 0; ell != numLambda; ++ell)
        {
            double lambda = lambdav[ell];
            for (int t = 0; t != numTheta; ++t)
            {
                double theta = thetav[t];
                S11vv(ell, t) = S11(lambda, theta);
                if (mode == ScatteringMode::SphericalPolarization || mode == ScatteringMode::SpheroidalPolarization)
                {
                    S12vv(ell, t) = S12(lambda, theta);
                    S33vv(ell, t) = S33(lambda, theta);
                    S34vv(ell, t) = S34(lambda, theta);
                }
            }
        }
    }

    return mu;
}

////////////////////////////////////////////////////////////////////

string SingleGrainDustMix::resourceNameForMuellerMatrix() const
{
    return string();
}

////////////////////////////////////////////////////////////////////
