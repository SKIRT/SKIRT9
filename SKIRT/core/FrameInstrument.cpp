/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FrameInstrument.hpp"
#include "FluxRecorder.hpp"
#include "PhotonPacket.hpp"

////////////////////////////////////////////////////////////////////

void FrameInstrument::setupSelfBefore()
{
    DistantInstrument::setupSelfBefore();

    // configure flux recorder
    instrumentFluxRecorder()->includeSurfaceBrightnessForDistant(
        numPixelsX(), numPixelsY(), fieldOfViewX() / numPixelsX(), fieldOfViewY() / numPixelsY(), centerX(), centerY());

    // precalculate information needed by pixelOnDetector() function
    _costheta = cos(inclination());
    _sintheta = sin(inclination());
    _cosphi = cos(azimuth());
    _sinphi = sin(azimuth());
    _cosomega = cos(roll());
    _sinomega = sin(roll());
    _Nxp = numPixelsX();
    _Nyp = numPixelsY();
    _xpmin = centerX() - 0.5 * fieldOfViewX();
    _xpsiz = fieldOfViewX() / numPixelsX();
    _ypmin = centerY() - 0.5 * fieldOfViewY();
    _ypsiz = fieldOfViewY() / numPixelsY();
}

////////////////////////////////////////////////////////////////////

void FrameInstrument::detect(PhotonPacket* pp)
{
    int l = pixelOnDetector(pp);
    instrumentFluxRecorder()->detect(pp, l);
}

////////////////////////////////////////////////////////////////////

int FrameInstrument::pixelOnDetector(const PhotonPacket* pp) const
{
    // get the position
    double x, y, z;
    pp->position().cartesian(x, y, z);

    // transform to detector coordinates using inclination, azimuth, and roll angle
    double xpp = -_sinphi * x + _cosphi * y;
    double ypp = -_cosphi * _costheta * x - _sinphi * _costheta * y + _sintheta * z;
    double xp = _cosomega * xpp - _sinomega * ypp;
    double yp = _sinomega * xpp + _cosomega * ypp;

    // scale and round to pixel index
    int i = static_cast<int>(floor((xp - _xpmin) / _xpsiz));
    int j = static_cast<int>(floor((yp - _ypmin) / _ypsiz));
    if (i < 0 || i >= _Nxp || j < 0 || j >= _Nyp)
        return -1;
    else
        return i + _Nxp * j;
}

////////////////////////////////////////////////////////////////////
