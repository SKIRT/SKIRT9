/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GeometricSource.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"
#include "WavelengthRangeInterface.hpp"

//////////////////////////////////////////////////////////////////////

int GeometricSource::dimension() const
{
    return _geometry->dimension();
}

//////////////////////////////////////////////////////////////////////

double GeometricSource::luminosity() const
{
    return _normalization->luminosity(_sed);
}

//////////////////////////////////////////////////////////////////////

double GeometricSource::specificLuminosity(double wavelength) const
{
    if (!interface<WavelengthRangeInterface>()->wavelengthRange().contains(wavelength)) return 0.;
    return _sed->specificLuminosity(wavelength) * luminosity();
}

//////////////////////////////////////////////////////////////////////

void GeometricSource::launch(PhotonPacket* pp, size_t historyIndex, double L) const
{
    // generate a random position from the geometry
    Position bfr = _geometry->generatePosition();

    // generate a random wavelength from the SED and/or from the bias distribution
    double lambda, w;
    double xi = wavelengthBias();
    if (!xi)
    {
        // no biasing -- simply use the intrinsic spectral distribution
        lambda = _sed->generateWavelength();
        w = 1.;
    }
    else
    {
        // biasing -- use one or the other distribution
        if (random()->uniform() > xi) lambda = sed()->generateWavelength();
        else lambda = wavelengthBiasDistribution()->generateWavelength();

        // calculate the compensating weight factor
        double s = sed()->specificLuminosity(lambda);
        if (!s)
        {
            // if the wavelength can't occur in the intrinsic distribution,
            // the weight factor is zero regardless of the probability in the bias distribution
            // (handling this separately also avoids NaNs in the pathological case s=d=0)
            w = 0.;
        }
        else
        {
            // regular composite bias weight
            double d = wavelengthBiasDistribution()->probability(lambda);
            w = s / ((1-xi)*s + xi*d);
        }
    }

    // launch the photon packet with isotropic direction
    pp->launch(historyIndex, lambda, L*w, bfr, random()->direction());
}

//////////////////////////////////////////////////////////////////////
