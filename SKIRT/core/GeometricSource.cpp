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
    double lambda = _sed->generateWavelength();
    // lambda = wavelengthBiasDistribution()->generateWavelength();

    // launch the photon packet with isotropic direction
    pp->launch(historyIndex, lambda, L, bfr, random()->direction());
}

//////////////////////////////////////////////////////////////////////
