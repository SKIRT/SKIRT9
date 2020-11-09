/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ReadFits3DGeometry.hpp"
#include "FITSInOut.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void ReadFits3DGeometry::setupSelfBefore()
{
    GenGeometry::setupSelfBefore();

    // Read in the datacube
    FITSInOut::read(this, _filename, _datacube, _nx, _ny, _nz);

    // Verify maximum size:
    //   we use NR::locate() on the cumulative distribution of the pixels, which returns an int,
    //   so things will fail badly if the datacube is larger than what fits in an int
    if (_datacube.size() > static_cast<size_t>(std::numeric_limits<int>::max() - 1))
        throw FATALERROR("Datacube size exceeds maximum supported size");

    // Construct a vector with the normalized cumulative density distribution
    double sum = NR::cdf(_Xv, _datacube);

    // Normalize the density datacube
    _datacube /= sum * _pixelScale * _pixelScale * _pixelScale;

    // Calculate the boundaries of the datacube in physical coordinates
    _xmin = -(_nx / 2.) * _pixelScale;
    _ymin = -(_ny / 2.) * _pixelScale;
    _zmin = -(_nz / 2.) * _pixelScale;
}

////////////////////////////////////////////////////////////////////

double ReadFits3DGeometry::density(Position bfr) const
{
    double x, y, z;
    bfr.cartesian(x, y, z);

    // Find the corresponding spaxel in the datacube
    int i = static_cast<int>(floor((x - _xmin) / _pixelScale));
    int j = static_cast<int>(floor((y - _ymin) / _pixelScale));
    int k = static_cast<int>(floor((z - _zmin) / _pixelScale));
    if (i < 0 || i >= _nx || j < 0 || j >= _ny || k < 0 || k >= _nz) return 0.;

    // Return the density
    return _datacube[i + _nx * j + _nx * _ny * k];
}

////////////////////////////////////////////////////////////////////

Position ReadFits3DGeometry::generatePosition() const
{
    // Draw a random position in cube based on the cumulative density distribution
    int index = NR::locate(_Xv, random()->uniform());
    int i = (index % (_nx * _ny)) % _nx;
    int j = ((index - i) % (_nx * _ny)) / _nx;
    int k = (index - _nx * j - i) / (_nx * _ny);

    // Determine the x, y and z coordinate in the datacube
    double x = _xmin + (i + random()->uniform()) * _pixelScale;
    double y = _ymin + (j + random()->uniform()) * _pixelScale;
    double z = _zmin + (k + random()->uniform()) * _pixelScale;

    // Return the position
    return Position(x, y, z);
}

////////////////////////////////////////////////////////////////////

double ReadFits3DGeometry::SigmaX() const
{
    double sum = 0;

    // For each position, get the density and add it to the total
    int j = static_cast<int>(floor(_ny / 2.));
    int k = static_cast<int>(floor(_nz / 2.));
    for (int i = 0; i < _nx; i++)
    {
        sum += _datacube[i + _nx * j + _nx * _ny * k];
    }

    // Return the x-axis surface density
    return sum * _pixelScale;
}

////////////////////////////////////////////////////////////////////

double ReadFits3DGeometry::SigmaY() const
{
    double sum = 0;

    // For each position, get the density and add it to the total
    int i = static_cast<int>(floor(_nx / 2.));
    int k = static_cast<int>(floor(_nz / 2.));
    for (int j = 0; j < _ny; j++)
    {
        sum += _datacube[i + _nx * j + _nx * _ny * k];
    }

    // Return the y-axis surface density
    return sum * _pixelScale;
}

////////////////////////////////////////////////////////////////////

double ReadFits3DGeometry::SigmaZ() const
{
    double sum = 0;

    // For each position, get the density and add it to the total
    int i = static_cast<int>(floor(_nx / 2.));
    int j = static_cast<int>(floor(_ny / 2.));
    for (int k = 0; k < _nz; k++)
    {
        sum += _datacube[i + _nx * j + _nx * _ny * k];
    }

    // Return the z-axis surface density
    return sum * _pixelScale;
}

////////////////////////////////////////////////////////////////////
