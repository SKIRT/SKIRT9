/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SphereSpatialGrid.hpp"
#include "FatalError.hpp"
#include "Mesh.hpp"

//////////////////////////////////////////////////////////////////////

void SphereSpatialGrid::setupSelfBefore()
{
    SpatialGrid::setupSelfBefore();
    if (_maxRadius <= _minRadius) throw FATALERROR("The outer radius must be larger than the inner radius");
}

//////////////////////////////////////////////////////////////////////

Box SphereSpatialGrid::boundingBox() const
{
    return Box(-_maxRadius, -_maxRadius, -_maxRadius, _maxRadius, _maxRadius, _maxRadius);
}

//////////////////////////////////////////////////////////////////////

int SphereSpatialGrid::initPolarGrid(Array& thetav, Array& cosv, const Mesh* mesh)
{
    // set up the polar grid; pre-calculate the cosines for each angular boundary
    double numTheta = mesh->numBins();
    thetav = mesh->mesh() * M_PI;
    cosv = cos(thetav);

    // make sure that the outer boundary cosine values are exact
    cosv[0] = 1.;
    cosv[numTheta] = -1.;

    // if there is a cosine value close to zero, make it exactly zero
    int numZeroes = 0;
    for (int j = 1; j < numTheta; j++)
    {
        if (abs(cosv[j]) < 1e-9)
        {
            numZeroes++;
            cosv[j] = 0.;
        }
    }
    if (numZeroes > 1) throw FATALERROR("There are multiple polar grid points very close to the equatorial plane");

    // if there is no cosine value close to zero, add an extra boundary
    if (numZeroes == 0)
    {
        // make a temporary copy of the original grid
        Array or_thetav = thetav;
        Array or_cv = cosv;

        // resize the grid to contain one extra bin; this clears all values
        numTheta++;
        thetav.resize(numTheta + 1);
        cosv.resize(numTheta + 1);

        // initialize the new grid with values of pi/2 so we can skip that index while copying
        thetav = M_PI_2;

        // copy the values from the original to the new grid, skipping the equatorial plane
        for (int j = 0; j < numTheta; j++)
        {
            int target = (or_cv[j] > 0) ? j : j + 1;
            thetav[target] = or_thetav[j];
            cosv[target] = or_cv[j];
        }
    }

    return numTheta;
}

//////////////////////////////////////////////////////////////////////
