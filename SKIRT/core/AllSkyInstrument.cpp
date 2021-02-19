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
    Vec co(_Ox - _Cx, _Oy - _Cy, _Oz - _Cz);  // vector in direction from crosshair to observer
    Vec up(_Ux, _Uy, _Uz);                    // vector in the upward direction
    if (co.norm() < 1e-20) throw FATALERROR("Crosshair is too close to observer");
    if (up.norm() < 1e-20) throw FATALERROR("Upwards direction cannot be null vector");
    if (Vec::cross(co, up).norm() < 1e-20)
        throw FATALERROR("Upwards direction cannot be parallel to viewing direction");

    // set number of pixels using fixed aspect ratio
    _Nx = 2 * _numPixelsY;
    _Ny = _numPixelsY;

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

    // rotate the axes into the alignment appropriate for our purposes (z-axis up, x-axis towards crosshair)
    _transform.rotateX(0., -1.);
    _transform.rotateZ(0., 1.);

    // determine the solid angle corresponding to each pixel, assuming an area preserving projection
    // and a usage fraction in the output rectangle of pi/4
    double omega = 16. / (_Nx * _Ny);

    // determine the scale of the output map axes, i.e. in geographical coordinates
    double inc = M_PI / _Ny;

    // configure flux recorder
    instrumentFluxRecorder()->includeSurfaceBrightnessForLocal(_Nx, _Ny, omega, -inc, inc, 0., 0., "posangle");
}

////////////////////////////////////////////////////////////////////

void AllSkyInstrument::determineSameObserverAsPreceding(const Instrument* precedingInstrument)
{
    auto other = dynamic_cast<const AllSkyInstrument*>(precedingInstrument);
    if (other && radius() == other->radius() && observerX() == other->observerX() && observerY() == other->observerY()
        && observerZ() == other->observerZ() && crossX() == other->crossX() && crossY() == other->crossY()
        && crossZ() == other->crossZ() && upX() == other->upX() && upY() == other->upY() && upZ() == other->upZ())
    {
        setSameObserverAsPreceding();
    }
}

////////////////////////////////////////////////////////////////////

Direction AllSkyInstrument::bfkobs(const Position& bfr) const
{
    // vector and distance from launch to observer
    Vec k = Vec(_Ox, _Oy, _Oz) - bfr;
    double d = k.norm();

    // if the distance is very small, return something arbitrary - the photon packet will not be detected anyway
    if (d <= _radius) return Direction();

    // otherwise return a unit vector in the direction from launch to observer
    return Direction(k / d);
}

////////////////////////////////////////////////////////////////////

Direction AllSkyInstrument::bfky(const Position& bfr) const
{
    // vector and distance from launch to observer
    Vec k = Vec(_Ox, _Oy, _Oz) - bfr;
    double d = k.norm();

    // if the distance is very small, return something arbitrary - the photon packet will not be detected anyway
    if (d <= _radius) return Direction();

    // vector in the plane normal to the line launch-observer
    // oriented along the projection of the up direction in that plane
    Vec ku(_Ux, _Uy, _Uz);
    Vec ky = Vec::cross(k, Vec::cross(ku, k));

    // return unit vector along y-axis
    return Direction(ky / ky.norm());
}

////////////////////////////////////////////////////////////////////

void AllSkyInstrument::detect(PhotonPacket* pp)
{
    // transform launch position from world to observer coordinates
    Position p(_transform.transform(pp->position()));

    // get the spherical coordinates of the launch position relative to the observer
    double d, inc, azi;
    p.spherical(d, inc, azi);

    // if the radial distance is very small, ignore the photon packet
    if (d <= _radius) return;

    // convert spherical coordinates to viewport coordinates:  -1 < x < 1  and -1 < y < 1
    double x, y;
    projection()->fromSphereToRectangle(inc, azi, x, y);

    // convert viewport coordinates to pixel indices
    int i = max(0, min(static_cast<int>((x + 1) * _Nx / 2.), _Nx - 1));
    int j = max(0, min(static_cast<int>((y + 1) * _Ny / 2.), _Ny - 1));

    // detect the photon packet in the appropriate pixel of the data cube and at the appropriate distance
    int l = i + _Nx * j;
    instrumentFluxRecorder()->detect(pp, l, d);
}

////////////////////////////////////////////////////////////////////
