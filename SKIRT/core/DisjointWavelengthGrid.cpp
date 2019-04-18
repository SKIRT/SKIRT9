/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DisjointWavelengthGrid.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

void DisjointWavelengthGrid::setupSelfAfter()
{
    WavelengthGrid::setupSelfAfter();

    // the subclass should have added wavelengths
    if (!_lambdav.size()) throw FATALERROR("Wavelength grid should have been initialized by subclass");
}

////////////////////////////////////////////////////////////////////

void DisjointWavelengthGrid::setWavelengthRange(const Array& lambdav, bool logScale)
{
    // copy and sort the specified characteristic wavelengths
    _lambdav = lambdav;
    std::sort(begin(_lambdav), end(_lambdav));
    size_t n = _lambdav.size();

    // verify that there is at least one wavelength and that the smallest one is positive
    if (!n) throw FATALERROR("There must be at least one wavelength in the grid");
    if (_lambdav[0] <= 0.0) throw FATALERROR("All wavelengths should be positive");

    // verify that there are no duplicates
    if (std::unique(begin(_lambdav), end(_lambdav)) != end(_lambdav))
        throw FATALERROR("There should be no duplicate wavelengths in the grid");

    // calculate the bin borders
    _lambdaleftv.resize(n);
    _lambdarightv.resize(n);
    _borderv.resize(n+1);
    if (n==1)
    {
        // just a single wavelength -> form narrow bin
        _lambdaleftv[0] = _borderv[0] = _lambdav[0] * 0.999;
        _lambdarightv[0] = _borderv[1] = _lambdav[0] * 1.001;
    }
    else
    {
        if (logScale)
        {
            _lambdaleftv[0] = _borderv[0] = sqrt(_lambdav[0]*_lambdav[0]*_lambdav[0]/_lambdav[1]);
            for (size_t ell=1; ell!=n; ++ell)
            {
                _lambdarightv[ell-1] = _lambdaleftv[ell] = _borderv[ell] = sqrt(_lambdav[ell-1]*_lambdav[ell]);
            }
            _lambdarightv[n-1] = _borderv[n] = sqrt(_lambdav[n-1]*_lambdav[n-1]*_lambdav[n-1]/_lambdav[n-2]);
        }
        else
        {
            _lambdaleftv[0] = _borderv[0] = (3.*_lambdav[0]-_lambdav[1])/2.;
            for (size_t ell=1; ell!=n; ++ell)
            {
                _lambdarightv[ell-1] = _lambdaleftv[ell] = _borderv[ell] = (_lambdav[ell-1]+_lambdav[ell])/2.;
            }
            _lambdarightv[n-1] = _borderv[n] = (3.*_lambdav[n-1]-_lambdav[n-2])/2.;
        }
    }

    // calculate the bin widths
    _dlambdav = _lambdarightv - _lambdaleftv;

    // setup the mapping from border bin indices to actual wavelength bin indices (see ell() function)
    _ellv.resize(n+2);
    _ellv[0] = -1;
    for (size_t ell=0; ell!=n; ++ell) _ellv[ell+1] = ell;
    _ellv[n+1] = -1;
}

////////////////////////////////////////////////////////////////////

void DisjointWavelengthGrid::setWavelengthBins(const Array& lambdav, double relativeHalfWidth, bool constantWidth)
{
    // copy and sort the specified characteristic wavelengths
    _lambdav = lambdav;
    std::sort(begin(_lambdav), end(_lambdav));
    size_t n = _lambdav.size();

    // verify that there is at least one wavelength and that the smallest one is positive
    if (!n) throw FATALERROR("There must be at least one wavelength in the grid");
    if (_lambdav[0] <= 0) throw FATALERROR("All wavelengths should be positive");

    // verify that the specified relative half width is positive
    if (relativeHalfWidth <= 0) throw FATALERROR("The relative wavelength bin width should be positive");

    // calculate the bin borders
    _lambdaleftv.resize(n);
    _lambdarightv.resize(n);
    _borderv.resize(2*n);
    if (!constantWidth)
    {
        for (size_t ell=0; ell!=n; ++ell)
        {
            _borderv[2*ell] = _lambdaleftv[ell] = _lambdav[ell] * (1.-relativeHalfWidth);
            _borderv[2*ell+1] = _lambdarightv[ell] = _lambdav[ell] * (1.+relativeHalfWidth);
        }
    }
    else
    {
        double delta = _lambdav[0] * relativeHalfWidth;
        for (size_t ell=0; ell!=n; ++ell)
        {
            _borderv[2*ell] = _lambdaleftv[ell] = _lambdav[ell] - delta;
            _borderv[2*ell+1] = _lambdarightv[ell] = _lambdav[ell] + delta;
        }
    }

    // verify that the bins do not overlap
    if (!std::is_sorted(begin(_borderv), end(_borderv)))
        throw FATALERROR("Non-adjacent wavelength bins should not overlap");

    // calculate the bin widths
    _dlambdav = _lambdarightv - _lambdaleftv;

    // setup the mapping from border bin indices to actual wavelength bin indices (see ell() function)
    _ellv.resize(2*n+1);
    _ellv[0] = -1;
    for (size_t ell=0; ell!=n; ++ell)
    {
        _ellv[2*ell+1] = ell;
        _ellv[2*ell+2] = -1;        // regions between the bins are considered out of range
    }
}

////////////////////////////////////////////////////////////////////

int DisjointWavelengthGrid::numBins() const
{
    return _lambdav.size();
}

////////////////////////////////////////////////////////////////////

double DisjointWavelengthGrid::wavelength(int ell) const
{
    return _lambdav[ell];
}

////////////////////////////////////////////////////////////////////

double DisjointWavelengthGrid::leftBorder(int ell) const
{
    return _lambdaleftv[ell];
}

////////////////////////////////////////////////////////////////////

double DisjointWavelengthGrid::rightBorder(int ell) const
{
    return _lambdarightv[ell];
}

////////////////////////////////////////////////////////////////////

double DisjointWavelengthGrid::effectiveWidth(int ell) const
{
    return _dlambdav[ell];
}

////////////////////////////////////////////////////////////////////

double DisjointWavelengthGrid::transmission(int /*ell*/, double /*lambda*/) const
{
    return 1.;
}

////////////////////////////////////////////////////////////////////

vector<int> DisjointWavelengthGrid::bins(double lambda) const
{
    // get the bin index, or -1 for "out of range"
    int ell = bin(lambda);

    // wrap the result in a list
    vector<int> result;
    if (ell >= 0) result.push_back(ell);
    return result;
}

////////////////////////////////////////////////////////////////////

int DisjointWavelengthGrid::bin(double lambda) const
{
    // get the index of the phantom wavelength bin defined by the list of all K borders (where K=N+1 or K=N*2)
    //  0  => out of range on the left side
    //  K  => out of range on the right side
    size_t index = std::upper_bound(begin(_borderv), end(_borderv), lambda) - begin(_borderv);

    // map this index to the actual wavelength bin index, or to -1 for "out of range"
    return _ellv[index];
}

////////////////////////////////////////////////////////////////////
