/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "NormalizedSource.hpp"
#include "Configuration.hpp"
#include "ContSED.hpp"
#include "FatalError.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void NormalizedSource::setupSelfBefore()
{
    Source::setupSelfBefore();

    auto config = find<Configuration>();
    _oligochromatic = config->oligochromatic();

    // cast the SED to a version that supports specific luminosities if possible
    _contsed = dynamic_cast<ContSED*>(_sed);

    // cache wavelength information depending on whether this is an oligo- or panchromatic simulation
    if (_oligochromatic)
    {
        _xi = 1.;  // always use bias distribution for oligochromatic simulations
        _biasDistribution = config->oligoWavelengthBiasDistribution();
    }
    else
    {
        _xi = _contsed ? wavelengthBias() : 0.;  // wavelength bias mechanism needs specific luminosities
        _biasDistribution = wavelengthBiasDistribution();
    }
}

//////////////////////////////////////////////////////////////////////

Range NormalizedSource::wavelengthRange() const
{
    return sed()->normalizationWavelengthRange();
}

//////////////////////////////////////////////////////////////////////

double NormalizedSource::luminosity() const
{
    return _normalization->luminosity(_sed);
}

//////////////////////////////////////////////////////////////////////

double NormalizedSource::specificLuminosity(double wavelength) const
{
    if (!_contsed) throw FATALERROR("Cannot probe specific luminosity for a line emission spectrum");

    if (!sed()->normalizationWavelengthRange().containsFuzzy(wavelength)) return 0.;
    return _contsed->specificLuminosity(wavelength) * luminosity();
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
        double s = _contsed->specificLuminosity(lambda);
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
    launchNormalized(pp, historyIndex, lambda, L * w);
}

//////////////////////////////////////////////////////////////////////
