/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Sphere1DSpatialGrid.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "PathSegmentGenerator.hpp"
#include "Random.hpp"
#include "SpatialGridPath.hpp"
#include "SpatialGridPlotFile.hpp"

//////////////////////////////////////////////////////////////////////

void Sphere1DSpatialGrid::setupSelfAfter()
{
    // Set up the grid properties
    _Nr = _meshRadial->numBins();
    _rv = minRadius() + _meshRadial->mesh() * (maxRadius() - minRadius());

    // base class setupSelfAfter() depends on initialization performed above
    SphereSpatialGrid::setupSelfAfter();
}

//////////////////////////////////////////////////////////////////////

int Sphere1DSpatialGrid::dimension() const
{
    return 1;
}

//////////////////////////////////////////////////////////////////////

int Sphere1DSpatialGrid::numCells() const
{
    return _Nr;
}

//////////////////////////////////////////////////////////////////////

double Sphere1DSpatialGrid::volume(int m) const
{
    int i = m;
    if (i < 0 || i >= _Nr)
        return 0.0;
    else
    {
        double rL = _rv[i];
        double rR = _rv[i + 1];
        return 4.0 * M_PI / 3.0 * (rR - rL) * (rR * rR + rR * rL + rL * rL);
    }
}

//////////////////////////////////////////////////////////////////////

int Sphere1DSpatialGrid::cellIndex(Position bfr) const
{
    return NR::locateFail(_rv, bfr.radius());
}

//////////////////////////////////////////////////////////////////////

Position Sphere1DSpatialGrid::centralPositionInCell(int m) const
{
    int i = m;
    double r = (_rv[i] + _rv[i + 1]) / 2.0;
    return Position(r, 0, 0);
}

//////////////////////////////////////////////////////////////////////

Position Sphere1DSpatialGrid::randomPositionInCell(int m) const
{
    int i = m;
    Direction bfk = random()->direction();
    double r = _rv[i] + (_rv[i + 1] - _rv[i]) * random()->uniform();
    return Position(r, bfk);
}

//////////////////////////////////////////////////////////////////////

class Sphere1DSpatialGrid::MySegmentGenerator : public PathSegmentGenerator
{
    const Sphere1DSpatialGrid* _grid{nullptr};
    enum class Phase { Inwards, Outwards };
    Phase _phase{Phase::Inwards};
    double _p{0.}, _q{0.};
    int _i{-1}, _imin{-1};

public:
    MySegmentGenerator(const Sphere1DSpatialGrid* grid) : _grid(grid) {}

    bool next() override
    {
        switch (state())
        {
            case State::Unknown:
            {
                // try moving the photon packet inside the grid
                double rmax = _grid->maxRadius();
                double rmin = _grid->minRadius();
                double r = sqrt(rx() * rx() + ry() * ry() + rz() * rz());
                _q = rx() * kx() + ry() * ky() + rz() * kz();
                _p = sqrt((r - _q) * (r + _q));
                if (r > rmax)
                {
                    if (_q > 0. || _p > rmax)
                    {
                        // there is no intersection with the grid; return an empty path
                        setState(State::Outside);
                        return false;
                    }
                    else
                    {
                        // path intersects rmax boundary going inward, qmax is negative; return the first empty segment
                        double qmax = -sqrt((rmax - _p) * (rmax + _p));
                        setEmptySegment(qmax - _q);
                        _i = _grid->_Nr - 1;
                        _q = qmax;
                        setState(State::Inside);
                        return true;
                    }
                }
                else if (r < rmin)
                {
                    // path intersects rmin boundary going outward, qmin is positive; return the first empty segment
                    double qmin = sqrt((rmin - _p) * (rmin + _p));
                    setEmptySegment(qmin - _q);
                    _i = 0;
                    _q = qmin;
                    setState(State::Inside);
                    return true;
                }

                // path starts inside the grid
                setState(State::Inside);

                // determine the initial grid cell
                _i = NR::locateClip(_grid->_rv, r);

                // determine the initial direction of movement relative to the center
                _phase = Phase::Outwards;
                if (_q < 0.)
                {
                    _imin = NR::locate(_grid->_rv, _p);  // returns -1 when p < rmin
                    if (_i > _imin) _phase = Phase::Inwards;
                }

                // fall through to determine the first actual segment
            }

            case State::Inside:
            {
                switch (_phase)
                {
                    case Phase::Inwards:
                    {
                        double rN = _grid->_rv[_i];  // i >= 0 here
                        double qN = -sqrt((rN - _p) * (rN + _p));
                        setSegment(_i, qN - _q);
                        _i--;  // i can become -1 here
                        _q = qN;
                        if (_i <= _imin) _phase = Phase::Outwards;
                        return true;
                    }
                    case Phase::Outwards:
                    {
                        double rN = _grid->_rv[_i + 1];
                        double qN = sqrt((rN - _p) * (rN + _p));
                        setSegment(_i, qN - _q);  // i can be -1 here, as intended
                        _i++;
                        _q = qN;
                        if (_i >= _grid->_Nr) setState(State::Outside);
                        return true;
                    }
                }
            }

            case State::Outside:
            {
            }
        }
        return false;
    }
};

/*
{
    // Determination of the initial position and direction of the path,
    // and calculation of some initial values

    path->clear();
    double x, y, z;
    path->position().cartesian(x, y, z);
    double kx, ky, kz;
    path->direction().cartesian(kx, ky, kz);
    double rmax = maxRadius();
    double rmin = minRadius();

    // Move the photon packet to the first grid cell that it will pass.
    // If it does not pass any grid cell, return an empty path.

    double r = path->position().radius();
    double q = x * kx + y * ky + z * kz;
    double p = sqrt((r - q) * (r + q));
    if (r > rmax)
    {
        if (q > 0.0 || p > rmax)
            return path->clear();
        else
        {
            r = rmax - 1e-8 * (_rv[_Nr] - _rv[_Nr - 1]);
            // path intersects rmax boundary going inward, qmax is negative
            double qmax = -sqrt((rmax - p) * (rmax + p));
            double ds = (qmax - q);
            path->addSegment(-1, ds);
            q = qmax;
        }
    }
    else if (r < rmin)
    {
        r = rmin + 1e-8 * (_rv[1] - _rv[0]);
        // path intersects rmin boundary going outward, qmin is positive
        double qmin = sqrt((rmin - p) * (rmin + p));
        double ds = (qmin - q);
        path->addSegment(-1, ds);
        q = qmin;
    }

    // Determination of the initial grid cell

    int i = NR::locateClip(_rv, r);

    // And here we go...

    double rN, qN;

    // Inward movement (not always...)

    if (q < 0.0)
    {
        int imin = NR::locate(_rv, p);  // returns -1 when p < rmin
        while (i > imin)
        {
            rN = _rv[i];  // i >= 0 here
            qN = -sqrt((rN - p) * (rN + p));
            int m = i;
            double ds = qN - q;
            path->addSegment(m, ds);
            i--;  // i can become -1 here
            q = qN;
        }
    }

    // Outward movement (starting with i potentially -1)

    while (true)
    {
        rN = _rv[i + 1];
        qN = sqrt((rN - p) * (rN + p));
        int m = i;
        double ds = qN - q;
        path->addSegment(m, ds);  // m can be -1 here, as intended
        i++;
        if (i >= _Nr)
            return;
        else
            q = qN;
    }
}
*/

//////////////////////////////////////////////////////////////////////

std::unique_ptr<PathSegmentGenerator> Sphere1DSpatialGrid::createPathSegmentGenerator() const
{
    return std::make_unique<MySegmentGenerator>(this);
}

//////////////////////////////////////////////////////////////////////

void Sphere1DSpatialGrid::write_xy(SpatialGridPlotFile* outfile) const
{
    for (int i = 0; i <= _Nr; i++) outfile->writeCircle(_rv[i]);
}

//////////////////////////////////////////////////////////////////////
