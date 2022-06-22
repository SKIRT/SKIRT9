/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ClumpyGeometryDecorator.hpp"
#include "NR.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void ClumpyGeometryDecorator::setupSelfAfter()
{
    GenGeometry::setupSelfAfter();

    // generate the random positions of the clumps
    _clumpv.resize(_numClumps);
    if (_seed) random()->push(_seed);
    for (int i = 0; i < _numClumps; i++) _clumpv[i] = _geometry->generatePosition();
    if (_seed) random()->pop();

    // sort the vector with the positions of the clumps in increasing x-coordinates
    NR::sort(_clumpv);
}

////////////////////////////////////////////////////////////////////

double ClumpyGeometryDecorator::density(Position bfr) const
{
    double rhosmooth = (1.0 - _clumpFraction) * _geometry->density(bfr);
    if (_cutoffClumps && !rhosmooth) return 0.0;  // don't allow clumps outside of smooth distribution

    double rhoclumpy = 0.0;
    double Mclump = _clumpFraction / _numClumps;  // total mass per clump
    int istart = max(0, NR::locate(_clumpv, bfr - Vec(_clumpRadius, 0, 0)));
    int iend = max(0, NR::locate(_clumpv, bfr + Vec(_clumpRadius, 0, 0)));
    for (int i = istart; i <= iend; i++)
    {
        double u = (bfr - _clumpv[i]).norm() / _clumpRadius;
        rhoclumpy += Mclump * _smoothingKernel->density(u) / pow(_clumpRadius, 3);
    }

    return rhosmooth + rhoclumpy;
}

////////////////////////////////////////////////////////////////////

Position ClumpyGeometryDecorator::generatePosition() const
{
    // loop until an appropriate position has been found
    while (true)
    {
        double X = random()->uniform();
        if (X > _clumpFraction)
            return _geometry->generatePosition();
        else
        {
            // random clump number based on X
            int i = min(static_cast<int>((X / _clumpFraction) * _numClumps), _numClumps - 1);
            double u = _smoothingKernel->generateRadius();
            Direction bfk(random()->direction());
            Position bfr(_clumpv[i] + u * _clumpRadius * bfk);

            // reject positions outside of smooth distribution
            if (!_cutoffClumps || _geometry->density(bfr)) return bfr;
        }
    }
}

////////////////////////////////////////////////////////////////////

double ClumpyGeometryDecorator::SigmaX() const
{
    return _geometry->SigmaX();
}

////////////////////////////////////////////////////////////////////

double ClumpyGeometryDecorator::SigmaY() const
{
    return _geometry->SigmaY();
}

////////////////////////////////////////////////////////////////////

double ClumpyGeometryDecorator::SigmaZ() const
{
    return _geometry->SigmaZ();
}

////////////////////////////////////////////////////////////////////
