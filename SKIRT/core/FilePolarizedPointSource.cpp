/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FilePolarizedPointSource.hpp"
#include "AngularDistribution.hpp"
#include "Configuration.hpp"
#include "ContSED.hpp"
#include "PhotonPacket.hpp"
#include "PolarizationProfile.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void FilePolarizedPointSource::setupSelfBefore()
{
    Source::setupSelfBefore();

    // read the input file
    // construct the SED averaged over all angles from the input file
    // _sed = ...
    // construct angular distribution object for each wavelength point
    // construct polarization profile object for each wavelength point
    // adjust for symmetry axis orientation in both sets of objects

    // cache wavelength sampling parameters depending on whether this is an oligo- or panchromatic simulation
    auto config = find<Configuration>();
    if (config->oligochromatic())
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
    if (hasVelocity()) _bvi = this;
}

////////////////////////////////////////////////////////////////////

void FilePolarizedPointSource::setupSelfAfter()
{
    Source::setupSelfAfter();

    // warn the user if this source's intrinsic wavelength range does not fully cover the configured wavelength range
    informAvailableWavelengthRange(_sed->intrinsicWavelengthRange(), _sed->type());
}

////////////////////////////////////////////////////////////////////

int FilePolarizedPointSource::dimension() const
{
    if (positionX() || positionY() || symmetryX() || symmetryY()) return 3;
    if (hasVelocity() && (velocityX() || velocityY())) return 3;
    return 2;
}

////////////////////////////////////////////////////////////////////

bool FilePolarizedPointSource::hasVelocity() const
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

////////////////////////////////////////////////////////////////////

Vec FilePolarizedPointSource::velocity() const
{
    return Vec(velocityX(), velocityY(), velocityZ());
}

////////////////////////////////////////////////////////////////////

Range FilePolarizedPointSource::wavelengthRange() const
{
    return _sed->normalizationWavelengthRange();
}

////////////////////////////////////////////////////////////////////

double FilePolarizedPointSource::luminosity() const
{
    return _normalization->luminosityForSED(_sed);
}

////////////////////////////////////////////////////////////////////

double FilePolarizedPointSource::specificLuminosity(double wavelength) const
{
    return _sed->normalizationWavelengthRange().containsFuzzy(wavelength)
               ? _sed->specificLuminosity(wavelength) * luminosity()
               : 0.;
}

////////////////////////////////////////////////////////////////////

void FilePolarizedPointSource::launch(PhotonPacket* pp, size_t historyIndex, double L) const
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
            lambda = _sed->generateWavelength();
        else
            lambda = _biasDistribution->generateWavelength();

        // calculate the compensating weight factor
        double s = _sed->specificLuminosity(lambda);
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

    // get the source position
    Position bfr(positionX(), positionY(), positionZ());

    // get or construct angular distribution object for this wavelength point
    AngularDistribution* angularDistribution = nullptr;

    // get or construct polarization profile object for this wavelength point
    PolarizationProfile* polarizationProfile = nullptr;

    // launch the photon packet with the proper wavelength, weight, position, direction,
    // bulk veloclity, angular distribution and polarization profile
    pp->launch(historyIndex, lambda, L * w, bfr, angularDistribution->generateDirection(), _bvi, angularDistribution,
               polarizationProfile);
}

////////////////////////////////////////////////////////////////////
