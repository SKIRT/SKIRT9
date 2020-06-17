/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "HEALPixSkyInstrument.hpp"
#include "FatalError.hpp"
#include "FluxRecorder.hpp"
#include "Log.hpp"
#include "PhotonPacket.hpp"

////////////////////////////////////////////////////////////////////

void HEALPixSkyInstrument::setupSelfBefore()
{
    Instrument::setupSelfBefore();

    // verify consistency of vector properties
    Vec co(_Ox - _Cx, _Oy - _Cy, _Oz - _Cz);  // vector in direction from crosshair to observer
    Vec up(_Ux, _Uy, _Uz);                    // vector in the upward direction
    if (co.norm() < 1e-20) throw FATALERROR("Crosshair is too close to observer");
    if (up.norm() < 1e-20) throw FATALERROR("Upwards direction cannot be null vector");
    if (Vec::cross(co, up).norm() < 1e-20)
        throw FATALERROR("Upwards direction cannot be parallel to viewing direction");

    // compute the side length of the subdivision within a single HEALPix base pixel
    _Nside = 1 << _order;

    // we map the 12 * _Nside^2 HEALPix pixels to an almost square image with 16 * _Nside * (_Nside-1) pixels
    _Nx = 4 * _Nside;
    _Ny = 4 * _Nside - 1;

    // determine linear size of a single pixel
    // each HEALPix pixel has the same spherical area pi/(3*_Nside^2)
    // we assume that the pixels are squares with this surface area
    // due to the calibration of the FrameInstrument, we need to multiply this with the radius of the instrument
    _s = sqrt(M_PI / (3 * _Nside * _Nside)) * _radius;

    // setup the transformation from world to observer coordinates

    // translate to observer position
    _transform.translate(-_Ox, -_Oy, -_Oz);

    // unit vector in direction from crosshair to observer
    Vec kn(_Ox - _Cx, _Oy - _Cy, _Oz - _Cz);
    kn /= kn.norm();
    double a = kn.x();
    double b = kn.y();
    double c = kn.z();

    // rotate from world to observer coordinates just as for the perspective transformation
    double v = sqrt(b * b + c * c);
    if (v > 0.3)
    {
        _transform.rotateX(c / v, -b / v);
        _transform.rotateY(v, -a);
        double k = (b * b + c * c) * _Ux - a * b * _Uy - a * c * _Uz;
        double l = c * _Uy - b * _Uz;
        double u = sqrt(k * k + l * l);
        _transform.rotateZ(l / u, -k / u);
    }
    else
    {
        v = sqrt(a * a + c * c);
        _transform.rotateY(c / v, -a / v);
        _transform.rotateX(v, -b);
        double k = c * _Ux - a * _Uz;
        double l = (a * a + c * c) * _Uy - a * b * _Ux - b * c * _Uz;
        double u = sqrt(k * k + l * l);
        _transform.rotateZ(l / u, -k / u);
    }
    // rather than flipping the z-axis as is done for the perspective transformation,
    // rotate the axes into the alignment appropriate for our purposes (z-axis up, x-axis towards crosshair)
    _transform.rotateX(0., 1.);
    _transform.rotateZ(0., -1.);

    // configure flux recorder with a large distance relative to the pixel size so that atan(s/2d) = s/2d
    // and the default calibration can be easily corrected when detecting each individual photon packet
    instrumentFluxRecorder()->setRestFrameDistance(_s * 1e8);
    instrumentFluxRecorder()->includeSurfaceBrightness(_Nx, _Ny, _s, _s, 0, 0, false);
}

////////////////////////////////////////////////////////////////////

void HEALPixSkyInstrument::determineSameObserverAsPreceding(const Instrument* precedingInstrument)
{
    auto other = dynamic_cast<const HEALPixSkyInstrument*>(precedingInstrument);
    if (other && radius() == other->radius() && observerX() == other->observerX() && observerY() == other->observerY()
        && observerZ() == other->observerZ() && crossX() == other->crossX() && crossY() == other->crossY()
        && crossZ() == other->crossZ() && upX() == other->upX() && upY() == other->upY() && upZ() == other->upZ())
    {
        setSameObserverAsPreceding();
    }
}

////////////////////////////////////////////////////////////////////

Direction HEALPixSkyInstrument::bfkobs(const Position& bfr) const
{
    // vector and distance from launch to observer
    Vec k = Vec(_Ox, _Oy, _Oz) - bfr;
    double d = k.norm();

    // if the distance is very small, return something arbitrary - the photon packet will not be detected anyway
    if (d < _s) return Direction();

    // otherwise return a unit vector in the direction from launch to observer
    return Direction(k / d);
}

////////////////////////////////////////////////////////////////////

Direction HEALPixSkyInstrument::bfkx(const Position& bfr) const
{
    // vector and distance from launch to observer
    Vec k = Vec(_Ox, _Oy, _Oz) - bfr;
    double d = k.norm();

    // if the distance is very small, return something arbitrary - the photon packet will not be detected anyway
    if (d < _s) return Direction();

    // vector in the plane normal to the line launch-observer
    // oriented perpendicular to the projection of the up direction in that plane
    Vec ku(_Ux, _Uy, _Uz);
    Vec kx = Vec::cross(ku, k);

    // return unit vector along y-axis
    return Direction(kx / kx.norm());
}

////////////////////////////////////////////////////////////////////

Direction HEALPixSkyInstrument::bfky(const Position& bfr) const
{
    // vector and distance from launch to observer
    Vec k = Vec(_Ox, _Oy, _Oz) - bfr;
    double d = k.norm();

    // if the distance is very small, return something arbitrary - the photon packet will not be detected anyway
    if (d < _s) return Direction();

    // vector in the plane normal to the line launch-observer
    // oriented along the projection of the up direction in that plane
    Vec ku(_Ux, _Uy, _Uz);
    Vec ky = Vec::cross(k, Vec::cross(ku, k));

    // return unit vector along y-axis
    return Direction(ky / ky.norm());
}

////////////////////////////////////////////////////////////////////

void HEALPixSkyInstrument::detect(PhotonPacket* pp)
{
    // transform launch position from world to observer coordinates
    Position p(_transform.transform(pp->position()));

    // get the spherical coordinates of the launch position relative to the observer
    double d, theta, phi;
    p.spherical(d, theta, phi);
    phi += M_PI;

    // if the radial distance is very small, ignore the photon packet
    if (d < _s) return;

    // the HEALPix mapping algorithm expects theta in [0, pi] and phi in [0, 2*pi]
    if (theta < 0.) theta += M_PI;
    if (phi < 0.) phi += 2. * M_PI;

    // the code below was mostly copied from healpix_base.cc in the official HEALPix repository
    // it corresponds to the loc2pix function using the RING scheme
    // unlike the loc2pix function, we do not compute a full RING pixel index, but restrict
    // ourselves to a computation of the ring index and pixel-in-ring index, which we respectively
    // use as vertical and horizontal pixel index
    double z = cos(theta);
    double za = abs(z);
    double tt = fmod(2. * phi / M_PI, 4.);

    // i: pixel-in-ring index (horizontal pixel index)
    // j: ring index (vertical pixel index)
    int i, j;
    // first figure out if we are in the equatorial or polar region
    if (za <= 2. / 3.)
    {
        // equatorial region: all rings have 4*_Nside pixels
        double t1 = _Nside * (0.5 + tt);
        double t2 = 0.75 * _Nside * z;
        int jp = static_cast<int>(floor(t1 - t2));
        int jm = static_cast<int>(floor(t1 + t2));

        // we first compute the ring index as the offset in the equatorial region (j in [1, 2*_Nside+1])
        j = _Nside + 1 + jp - jm;
        int kshift = 1 - (j & 1);
        int temp = jp + jm + kshift + 1 + 7 * _Nside;
        i = (_order > 0) ? ((temp >> 1) & (4 * _Nside - 1)) : ((temp >> 1) % (4 * _Nside));
        // we now add the correct offset to the ring index to account for the rings in the north polar region
        j += _Nside - 2;
    }
    else
    {
        // polar region: the number of pixels per ring depends on how far the ring is from the pole
        double tp = tt - static_cast<int>(tt);
        double tmp = (za < 0.99) ? (_Nside * sqrt(3. * (1. - za))) : (_Nside * sin(theta) / sqrt((1. + za) / 3.));

        int jp = static_cast<int>(tp * tmp);
        int jm = static_cast<int>((1. - tp) * tmp);

        // we first compute the ring index as an offset from the pole (j in [1, _Nside-1])
        j = jp + jm + 1;
        i = static_cast<int>(tt * j);
        // we now convert the ring index to the actual ring index, and distinguish between north and south pole
        if (z < 0)
        {
            j = 4 * _Nside - j - 1;
        }
        else
        {
            --j;
        }
    }

    // detect the photon packet in the appropriate pixel of the data cube and at the appropriate distance
    int l = i + _Nx * j;
    instrumentFluxRecorder()->detect(pp, l, d);
}

////////////////////////////////////////////////////////////////////
