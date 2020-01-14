/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "NormalizedSource.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void NormalizedSource::setupSelfBefore()
{
    Source::setupSelfBefore();

    auto config = find<Configuration>();
    _oligochromatic = config->oligochromatic();

    // cache wavelength information depending on whether this is an oligo- or panchromatic simulation
    if (_oligochromatic)
    {
        _xi = 1.;  // always use bias distribution for oligochromatic simulations
        _biasDistribution = config->oligoWavelengthBiasDistribution();
    }
    else
    {
        _xi = wavelengthBias();
        _biasDistribution = wavelengthBiasDistribution();
    }

    // if we have a nonzero bulk velocity, set the interface object to ourselves; otherwise leave it at null pointer
    if (!_oligochromatic && (velocityX() || velocityY() || velocityZ())) _bvi = this;
}

//////////////////////////////////////////////////////////////////////

int NormalizedSource::dimension() const
{
    int velocityDimension = 1;
    if (!find<Configuration>()->oligochromatic())
    {
        if (velocityZ()) velocityDimension = 2;
        if (velocityX() || velocityY()) velocityDimension = 3;
    }
    return max(geometryDimension(), velocityDimension);
}

//////////////////////////////////////////////////////////////////////

bool NormalizedSource::hasVelocity() const
{
    if (velocityX() || velocityY() || velocityZ())
    {
        // refuse velocity for oligochromatic simulations
        // (this function is called from Configure so we cannot precompute this during setup)
        auto config = find<Configuration>();
        if (!config->oligochromatic()) return true;
    }
    return false;
}

//////////////////////////////////////////////////////////////////////

Range NormalizedSource::wavelengthRange() const
{
    return sed()->wavelengthRange();
}

//////////////////////////////////////////////////////////////////////

double NormalizedSource::luminosity() const
{
    return _normalization->luminosity(_sed);
}

//////////////////////////////////////////////////////////////////////

double NormalizedSource::specificLuminosity(double wavelength) const
{
    if (!sed()->wavelengthRange().containsFuzzy(wavelength)) return 0.;
    return _sed->specificLuminosity(wavelength) * luminosity();
}

//////////////////////////////////////////////////////////////////////

Vec NormalizedSource::velocity() const
{
    return Vec(velocityX(), velocityY(), velocityZ());
}

//////////////////////////////////////////////////////////////////////

void NormalizedSource::launch(PhotonPacket* pp, size_t historyIndex, double L) const
{
    // generate a random wavelength from the SED and/or from the bias distribution
    double lambda, w;
    if (!_xi)
    {
        // no biasing -- simply use the intrinsic spectral distribution
        lambda = _sed->generateWavelength();
        w = 1.;
    }
    else
    {
        // biasing -- use one or the other distribution
        if (random()->uniform() > _xi)
            lambda = sed()->generateWavelength();
        else
            lambda = _biasDistribution->generateWavelength();

        // calculate the compensating weight factor
        double s = sed()->specificLuminosity(lambda);
        if (!s)
        {
            // if the wavelength can't occur in the intrinsic distribution,
            // the weight factor is zero regardless of the probability in the bias distribution
            // (handling this separately also avoids NaNs in the pathological case s=b=0)
            w = 0.;
        }
        else
        {
            // regular composite bias weight
            double b = _biasDistribution->probability(lambda);
            w = s / ((1 - _xi) * s + _xi * b);
        }
    }

    // cause the subclass to launch the photon packet
    launchNormalized(pp, historyIndex, lambda, L * w, _bvi);
}

//////////////////////////////////////////////////////////////////////
