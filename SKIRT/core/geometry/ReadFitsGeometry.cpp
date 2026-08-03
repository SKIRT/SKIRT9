/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ReadFitsGeometry.hpp"
#include "FITSInOut.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void ReadFitsGeometry::setupSelfBefore()
{
    GenGeometry::setupSelfBefore();

    // Read the input file
    FITSInOut::read(this, _filename, _image, _nx, _ny, _nz);

    // Verify the image properties
    if (_numPixelsX != _nx || _numPixelsY != _ny)
        throw FATALERROR("Specified image pixel size does not correspond with the FITS image frame size");
    if (_nz != 1) throw FATALERROR("FITS image contains multiple frames");

    // Normalize the luminosities
    _image /= _image.sum();

    // Construct a vector with the normalized cumulative luminosities
    NR::cdf(_Xv, _image);

    // Calculate the boundaries of the image in physical coordinates
    _xmax = ((_numPixelsX - _centerX) * _pixelScale);
    _xmin = -_centerX * _pixelScale;
    _ymax = ((_numPixelsY - _centerY) * _pixelScale);
    _ymin = -_centerY * _pixelScale;

    // Calculate the sines and cosines of the position angle and inclination
    _cospa = cos(_positionAngle);
    _sinpa = sin(_positionAngle);
    _cosi = cos(_inclination);
    _sini = sin(_inclination);

    // Calculate the physical pixel size in the x direction of the galactic plane
    _deltax = _pixelScale / _cosi;

    // Calculate the coordinates of the 4 corners of the image in the rotated plane
    _C1x = _xmax;
    _C1y = _ymax;
    _C2x = _xmin;
    _C2y = _ymax;
    _C3x = _xmin;
    _C3y = _ymin;
    _C4x = _xmax;
    _C4y = _ymin;
    derotate(_C1x, _C1y);
    derotate(_C2x, _C2y);
    derotate(_C3x, _C3y);
    derotate(_C4x, _C4y);
}

////////////////////////////////////////////////////////////////////

double ReadFitsGeometry::density(Position bfr) const
{
    double x, y, z;
    bfr.cartesian(x, y, z);

    // Project and rotate the x and y coordinates
    project(x);
    rotate(x, y);

    // Find the corresponding pixel in the image
    int i = static_cast<int>(floor(x - _xmin) / _pixelScale);
    int j = static_cast<int>(floor(y - _ymin) / _pixelScale);
    if (i < 0 || i >= _numPixelsX || j < 0 || j >= _numPixelsY) return 0.0;

    // Return the density
    return _image[i + _nx * j] * exp(-fabs(z) / _scaleHeight) / (2 * _scaleHeight) / (_deltax * _pixelScale);
}

////////////////////////////////////////////////////////////////////

Position ReadFitsGeometry::generatePosition() const
{
    // Draw a random position in the plane of the galaxy based on the
    // cumulative luminosities per pixel
    double X1 = random()->uniform();
    int k = NR::locate(_Xv, X1);
    int i = k % _numPixelsX;
    int j = (k - i) / _numPixelsX;

    // Determine the x and y coordinate in the plane of the image
    double x = _xmin + (i + random()->uniform()) * _pixelScale;
    double y = _ymin + (j + random()->uniform()) * _pixelScale;

    // Derotate and deproject the x and y coordinates
    derotate(x, y);
    deproject(x);

    // Draw a random position along the minor axis
    double X2 = random()->uniform();
    double z = (X2 <= 0.5) ? _scaleHeight * log(2.0 * X2) : -_scaleHeight * log(2.0 * (1.0 - X2));

    // Return the position
    return Position(x, y, z);
}

////////////////////////////////////////////////////////////////////

double ReadFitsGeometry::SigmaX() const
{
    const int NSAMPLES = 10000;
    double sum = 0;

    // Find the maximum and minimum possible x value
    double xmax = max({_C1x, _C2x, _C3x, _C4x});
    double xmin = min({_C1x, _C2x, _C3x, _C4x});
    deproject(xmax);
    deproject(xmin);

    // For each position, get the density and add it to the total
    for (int k = 0; k < NSAMPLES; k++)
    {
        sum += density(Position(xmin + k * (xmax - xmin) / NSAMPLES, 0, 0));
    }

    // Return the x-axis surface density
    return (sum / NSAMPLES) * (xmax - xmin);
}

////////////////////////////////////////////////////////////////////

double ReadFitsGeometry::SigmaY() const
{
    const int NSAMPLES = 10000;
    double sum = 0;

    // Find the maximum and minimum possible y value
    double ymax = max({_C1y, _C2y, _C3y, _C4y});
    double ymin = min({_C1y, _C2y, _C3y, _C4y});

    // For each position, get the density and add it to the total
    for (int k = 0; k < NSAMPLES; k++)
    {
        sum += density(Position(0, ymin + k * (ymax - ymin) / NSAMPLES, 0));
    }

    // Return the y-axis surface density
    return (sum / NSAMPLES) * (ymax - ymin);
}

////////////////////////////////////////////////////////////////////

double ReadFitsGeometry::SigmaZ() const
{
    // Get the index of the luminosity vector for the center of the galaxy (-_xpmin,-_ypmin)
    int i = static_cast<int>(floor((-_xmin) / _pixelScale));
    int j = static_cast<int>(floor((-_ymin) / _pixelScale));

    // Return the z-axis surface density
    return _image[i + _nx * j] / (_pixelScale * _pixelScale);
}

////////////////////////////////////////////////////////////////////

void ReadFitsGeometry::rotate(double& x, double& y) const
{
    // Cache the original values of x and y
    double xorig = x;
    double yorig = y;

    // Calculate the coordinates in the plane of the image
    x = (_sinpa * xorig) + (_cospa * yorig);
    y = (-_cospa * xorig) + (_sinpa * yorig);
}

////////////////////////////////////////////////////////////////////

void ReadFitsGeometry::derotate(double& x, double& y) const
{
    // Cache the original values of x and y
    double xorig = x;
    double yorig = y;

    // Calculate the coordinates in the rotated plane
    x = (_sinpa * xorig) - (_cospa * yorig);
    y = (_cospa * xorig) + (_sinpa * yorig);
}

////////////////////////////////////////////////////////////////////

void ReadFitsGeometry::project(double& x) const
{
    x = x * _cosi;
}

////////////////////////////////////////////////////////////////////

void ReadFitsGeometry::deproject(double& x) const
{
    x = x / _cosi;
}

////////////////////////////////////////////////////////////////////
