/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ApertureInstrument.hpp"
#include "PhotonPacket.hpp"

////////////////////////////////////////////////////////////////////

void ApertureInstrument::setupSelfBefore()
{
    DistantInstrument::setupSelfBefore();

    // precalculate information needed by isInsideAperture() function
    _radius2 = radius() * radius();
    _costheta = cos(inclination());
    _sintheta = sin(inclination());
    _cosphi = cos(azimuth());
    _sinphi = sin(azimuth());
}

////////////////////////////////////////////////////////////////////

bool ApertureInstrument::isInsideAperture(PhotonPacket* pp) const
{
    // if the instrument has an aperture
    if (_radius2)
    {
        // get the photon packet position
        double x, y, z;
        pp->position().cartesian(x, y, z);

        // transform to detector coordinates using inclination and azimuth
        // but without performing the roll, which would not alter the radius
        double xpp = -_sinphi * x + _cosphi * y;
        double ypp = -_cosphi * _costheta * x - _sinphi * _costheta * y + _sintheta * z;

        // if the photon packet projects outside of the aperture, ignore it
        double radius2 = xpp * xpp + ypp * ypp;
        if (radius2 > _radius2) return false;
    }
    return true;
}

////////////////////////////////////////////////////////////////////
