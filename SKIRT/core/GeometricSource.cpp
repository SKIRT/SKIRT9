/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GeometricSource.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void GeometricSource::setupSelfBefore()
{
    Source::setupSelfBefore();

    _sed->setWavelengthRange(minWavelength(), maxWavelength());
}

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
    if (wavelength < minWavelength() || wavelength > maxWavelength()) return 0.;
    return _sed->specificLuminosity(wavelength) * luminosity();
}

//////////////////////////////////////////////////////////////////////

void GeometricSource::launch(PhotonPacket* pp, size_t historyIndex, double L) const
{
    // generate a random position from the geometry
    Position bfr = _geometry->generatePosition();

    // generate a random wavelength from the SED
    double lambda = _sed->generateWavelength();

    // launch the photon packet with isotropic direction
    pp->launch(historyIndex, lambda, L, bfr, random()->direction());
}

//////////////////////////////////////////////////////////////////////
