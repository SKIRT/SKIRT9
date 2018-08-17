/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AllSkyInstrument.hpp"
#include "FatalError.hpp"
#include "FluxRecorder.hpp"
#include "PhotonPacket.hpp"

////////////////////////////////////////////////////////////////////

void AllSkyInstrument::setupSelfBefore()
{
    Instrument::setupSelfBefore();

    // verify consistency of vector properties
    Vec co(_Ox-_Cx, _Oy-_Cy, _Oz-_Cz); // vector in direction from crosshair to observer
    Vec up(_Ux, _Uy, _Uz);             // vector in the upward direction
    if (co.norm() < 1e-20) throw FATALERROR("Crosshair is too close to observer");
    if (up.norm() < 1e-20) throw FATALERROR("Upwards direction cannot be null vector");
    if (Vec::cross(co,up).norm() < 1e-20) throw FATALERROR("Upwards direction cannot be parallel to viewing direction");

    // set number of pixels using fixed aspect ratio
    _Nx = 2 * _numPixelsY;
    _Ny = _numPixelsY;

    // determine linear size of a single pixel
    // assume: square pixels; fraction of effective pixels pi/4; total area sphere with radius R
    _s = sqrt(8)*_radius/_numPixelsY;

    // setup the transformation from world to observer coordinates

    // translate to observer position
    _transform.translate(-_Ox, -_Oy, -_Oz);

    // unit vector in direction from crosshair to observer
    Vec kn(_Ox-_Cx, _Oy-_Cy, _Oz-_Cz);
    kn /= kn.norm();
    double a = kn.x();
    double b = kn.y();
    double c = kn.z();

    // rotate from world to observer coordinates just as for the perspective transformation
    double v = sqrt(b*b+c*c);
    if (v > 0.3)
    {
        _transform.rotateX(c/v, -b/v);
        _transform.rotateY(v, -a);
        double k = (b*b+c*c)*_Ux - a*b*_Uy - a*c*_Uz;
        double l = c*_Uy - b*_Uz;
        double u = sqrt(k*k+l*l);
        _transform.rotateZ(l/u, -k/u);
    }
    else
    {
        v = sqrt(a*a+c*c);
        _transform.rotateY(c/v, -a/v);
        _transform.rotateX(v, -b);
        double k = c*_Ux - a*_Uz;
        double l = (a*a+c*c)*_Uy - a*b*_Ux - b*c*_Uz;
        double u = sqrt(k*k+l*l);
        _transform.rotateZ(l/u, -k/u);
    }
    // rather than flipping the z-axis as is done for the perspective transformation,
    // rotate the axes into the alignment appropriate for our purposes (z-axis up, x-axis towards crosshair)
    _transform.rotateX(0., 1.);
    _transform.rotateZ(0., -1.);

    // configure flux recorder with a large distance relative to the pixel size so that atan(s/2d) = s/2d
    // and the default calibration can be easily corrected when detecting each individual phoon packet
    double distance = _s * 1e8;
    instrumentFluxRecorder()->includeSurfaceBrightness(distance, _Nx, _Ny, _s, _s, 0, 0);
}

////////////////////////////////////////////////////////////////////

Direction AllSkyInstrument::bfkobs(const Position& bfr) const
{
    // vector and distance from launch to observer
    Vec k = Vec(_Ox,_Oy,_Oz) - bfr;
    double d = k.norm();

    // if the distance is very small, return something arbitrary - the photon packet will not be detected anyway
    if (d < _s) return Direction();

    // otherwise return a unit vector in the direction from launch to observer
    return Direction(k/d);
}

////////////////////////////////////////////////////////////////////

Direction AllSkyInstrument::bfkx(const Position& bfr) const
{
    // vector and distance from launch to observer
    Vec k = Vec(_Ox,_Oy,_Oz) - bfr;
    double d = k.norm();

    // if the distance is very small, return something arbitrary - the photon packet will not be detected anyway
    if (d < _s) return Direction();

    // vector in the plane normal to the line launch-observer
    // oriented perpendicular to the projection of the up direction in that plane
    Vec ku(_Ux, _Uy, _Uz);
    Vec kx = Vec::cross(ku,k);

    // return unit vector along y-axis
    return Direction(kx/kx.norm());
}

////////////////////////////////////////////////////////////////////

Direction AllSkyInstrument::bfky(const Position& bfr) const
{
    // vector and distance from launch to observer
    Vec k = Vec(_Ox,_Oy,_Oz) - bfr;
    double d = k.norm();

    // if the distance is very small, return something arbitrary - the photon packet will not be detected anyway
    if (d < _s) return Direction();

    // vector in the plane normal to the line launch-observer
    // oriented along the projection of the up direction in that plane
    Vec ku(_Ux, _Uy, _Uz);
    Vec ky = Vec::cross(k,Vec::cross(ku,k));

    // return unit vector along y-axis
    return Direction(ky/ky.norm());
}

////////////////////////////////////////////////////////////////////

void AllSkyInstrument::detect(PhotonPacket* pp)
{
    // transform launch position from world to observer coordinates
    Position p( _transform.transform(pp->position()) );

    // get the spherical coordinates of the launch position relative to the observer
    double d, inc, azi;
    p.spherical(d, inc, azi);

    // if the radial distance is very small, ignore the photon packet
    if (d < _s) return;

    // convert spherical coordinates to viewport coordinates:  -1 < x < 1  and -1 < y < 1
    double x, y;
    projection()->fromSphereToRectangle(inc, azi, x, y);

    // convert viewport coordinates to pixel indices
    int i = max(0, min(static_cast<int>((x+1)*_Nx/2.), _Nx-1));
    int j = max(0, min(static_cast<int>((y+1)*_Ny/2.), _Ny-1));

    // detect the photon packet in the appropriate pixel of the data cube and at the appropriate distance
    int l = i + _Nx*j;
    instrumentFluxRecorder()->detect(pp, l, d);
}

////////////////////////////////////////////////////////////////////
