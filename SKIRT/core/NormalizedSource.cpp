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

    if (_oligochromatic)
    {
        // cache the oligochromatic wavelengths
        _oligoWavelengths = config->oligoWavelengths();
        _oligoTotalBinWidth = config->oligoBinWidth() * _oligoWavelengths.size();
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

Vec NormalizedSource::bulkVelocity() const
{
    return Vec(velocityX(), velocityY(), velocityZ());
}

//////////////////////////////////////////////////////////////////////

void NormalizedSource::launch(PhotonPacket* pp, size_t historyIndex, double L) const
{
    double lambda, w;

    if (_oligochromatic)
    {
        // generate one of the oligochromatic wavelengths with equal probability
        size_t index = static_cast<size_t>(random()->uniform() * _oligoWavelengths.size());
        lambda = _oligoWavelengths[index];

        // bias the weight for the actual luminosity at this wavelength
        w = _oligoTotalBinWidth * sed()->specificLuminosity(lambda);
    }
    else
    {
        // generate a random wavelength from the SED and/or from the bias distribution
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
                // (handling this separately also avoids NaNs in the pathological case s=b=0)
                w = 0.;
            }
            else
            {
                // regular composite bias weight
                double b = wavelengthBiasDistribution()->probability(lambda);
                w = s / ((1-xi)*s + xi*b);
            }
        }
    }

    // cause the subclass to launch the photon packet
    launchNormalized(pp, historyIndex, lambda, L*w, _bvi);
}

//////////////////////////////////////////////////////////////////////
