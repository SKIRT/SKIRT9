/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ReadFits3DGeometry.hpp"
#include "FatalError.hpp"
#include "FITSInOut.hpp"
#include "NR.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void ReadFits3DGeometry::setupSelfBefore()
{
    GenGeometry::setupSelfBefore();

    // Read in the datacube
    FITSInOut::read(this, _filename, _datacube, _nx, _ny, _nz);

    // Normalize the datacube
    _datacube /= _datacube.sum();

    // Construct a vector with the normalized cumulative distribution
    NR::cdf(_Xv, _datacube);

    // Calculate the boundaries of the datacube in physical coordinates
    _xmin = -(_nx/2.)*_pixelScale;
    _xmax =  (_nx/2.)*_pixelScale;
    _ymin = -(_ny/2.)*_pixelScale;
    _ymax =  (_ny/2.)*_pixelScale;
    _zmin = -(_nz/2.)*_pixelScale;
    _zmax =  (_nz/2.)*_pixelScale;

}

////////////////////////////////////////////////////////////////////

double ReadFits3DGeometry::density(Position bfr) const
{
    double x,y,z;
    bfr.cartesian(x,y,z);

    // Find the corresponding spaxel in the datacube
    int i = static_cast<int>(floor((x-_xmin)/_pixelScale));
    int j = static_cast<int>(floor((y-_ymin)/_pixelScale));
    int k = static_cast<int>(floor((z-_zmin)/_pixelScale));
    if (i<0 || i>=_nx || j<0 || j>=_ny|| k<0 || k>=_nz) return 0.0;

    // Return the density
    return _datacube[i + _nx*j + _nx*_ny*k];
}

////////////////////////////////////////////////////////////////////

Position ReadFits3DGeometry::generatePosition() const
{
    // Draw a random position in cube based on the
    // cumulative distribution in each spaxel
    double X1 = random()->uniform();
    int Xindex = NR::locate(_Xv,X1);

    int i = (Xindex%(_nx*_ny))%_nx;
    int j = ((Xindex-i)%(_nx*_ny))/_nx;
    int k = (Xindex-_nx*j - i)/(_nx*_ny);

    // Determine the x, y and z coordinate in the datacube
    double x = _xmin + (i+random()->uniform())*_pixelScale;
    double y = _ymin + (j+random()->uniform())*_pixelScale;
    double z = _zmin + (k+random()->uniform())*_pixelScale;


    // Return the position
    return Position(x,y,z);
}

////////////////////////////////////////////////////////////////////

double ReadFits3DGeometry::SigmaX() const
{
    double sum = 0;
    const int NSAMPLES = 10000;
    double step = (_xmax-_xmin)/NSAMPLES;

    // For each position, get the density and add it to the total
    for (int k = 0; k < NSAMPLES; k++)
    {
        sum += density(Position(_xmin + k*step, 0, 0));
    }

    // Return the x-axis surface density
    return sum*step;
}

////////////////////////////////////////////////////////////////////

double ReadFits3DGeometry::SigmaY() const
{
    double sum = 0;
    const int NSAMPLES = 10000;
    double step = (_ymax-_ymin)/NSAMPLES;

    // For each position, get the density and add it to the total
    for (int k = 0; k < NSAMPLES; k++)
    {
        sum += density(Position(_ymin + k*step, 0, 0));
    }

    // Return the y-axis surface density
    return sum*step;
}

////////////////////////////////////////////////////////////////////

double ReadFits3DGeometry::SigmaZ() const
{
    double sum = 0;
    const int NSAMPLES = 10000;
    double step = (_zmax-_zmin)/NSAMPLES;

    // For each position, get the density and add it to the total
    for (int k = 0; k < NSAMPLES; k++)
    {
        sum += density(Position(0, _zmin + k*step, 0));
    }

    // Return the z-axis surface density
    return sum*step;
}

////////////////////////////////////////////////////////////////////
