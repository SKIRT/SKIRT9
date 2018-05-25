/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DistantInstrument.hpp"

////////////////////////////////////////////////////////////////////

void DistantInstrument::setupSelfBefore()
{
    Instrument::setupSelfBefore();

    // calculate sine and cosine for our angles
    double costheta = cos(_inclination);
    double sintheta = sin(_inclination);
    double cosphi = cos(_azimuth);
    double sinphi = sin(_azimuth);
    double cosomega = cos(_rollAngle);
    double sinomega = sin(_rollAngle);

    // calculate relevant directions
    _bfkobs = Direction(_inclination,_azimuth);
    _bfkx = Direction( + cosphi*costheta*sinomega - sinphi*cosomega,
                       + sinphi*costheta*sinomega + cosphi*cosomega,
                       - sintheta*sinomega );
    _bfky = Direction( - cosphi*costheta*cosomega - sinphi*sinomega,
                       - sinphi*costheta*cosomega + cosphi*sinomega,
                       + sintheta*cosomega );
}

////////////////////////////////////////////////////////////////////

Direction DistantInstrument::bfkobs(const Position& /*bfr*/) const
{
    return _bfkobs;
}

////////////////////////////////////////////////////////////////////

Direction DistantInstrument::bfkx(const Position& /*bfr*/) const
{
    return _bfkx;
}

////////////////////////////////////////////////////////////////////

Direction DistantInstrument::bfky(const Position& /*bfr*/) const
{
    return _bfky;
}

////////////////////////////////////////////////////////////////////
