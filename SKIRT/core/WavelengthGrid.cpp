/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "WavelengthGrid.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

void WavelengthGrid::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    // the subclass should have added wavelengths
    if (!_lambdav.size()) throw FATALERROR("Wavelength grid should have been initialized by subclass");
}

////////////////////////////////////////////////////////////////////

void WavelengthGrid::setWavelengthRange(const Array& lambdav)
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

    _lambdaleftv[0] = _borderv[0] = _lambdav[0] * 0.999;
    for (size_t ell=1; ell!=n; ++ell)
    {
        _lambdarightv[ell-1] = _lambdaleftv[ell] = _borderv[ell] = sqrt(_lambdav[ell-1]*_lambdav[ell]);
    }
    _lambdarightv[n-1] = _borderv[n] = _lambdav[n-1] * 1.001;

    // calculate the bin widths
    _dlambdav = _lambdarightv - _lambdaleftv;

    // setup the mapping from border bin indices to actual wavelength bin indices (see ell() function)
    _ellv.resize(n+2);
    _ellv[0] = -1;
    for (size_t ell=0; ell!=n; ++ell) _ellv[ell+1] = ell;
    _ellv[n+1] = -1;
}

////////////////////////////////////////////////////////////////////

void WavelengthGrid::setWavelengthBins(const Array& lambdav, double relativeHalfWidth)
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
    for (size_t ell=0; ell!=n; ++ell)
    {
        _borderv[2*ell] = _lambdaleftv[ell] = _lambdav[ell] * (1.-relativeHalfWidth);
        _borderv[2*ell+1] = _lambdarightv[ell] = _lambdav[ell] * (1.+relativeHalfWidth);
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

int WavelengthGrid::numWavelengths() const
{
    return _lambdav.size();
}

////////////////////////////////////////////////////////////////////

double WavelengthGrid::lambda(int ell) const
{
    return _lambdav[ell];
}

////////////////////////////////////////////////////////////////////

double WavelengthGrid::dlambda(int ell) const
{
    return _dlambdav[ell];
}

////////////////////////////////////////////////////////////////////

double WavelengthGrid::lambdaLeft(int ell) const
{
    return _lambdaleftv[ell];
}

////////////////////////////////////////////////////////////////////

double WavelengthGrid::lambdaRight(int ell) const
{
    return _lambdarightv[ell];
}

////////////////////////////////////////////////////////////////////

int WavelengthGrid::ell(double lambda) const
{
    // get the index of the phantom wavelength bin defined by the list of all K borders (where K=N+1 or K=N*2)
    //  0  => out of range on the left side
    //  K  => out of range on the right side
    size_t index = std::upper_bound(begin(_borderv), end(_borderv), lambda) - begin(_borderv);

    // map this index to the actual wavelength bin index, or to "out of range"
    return _ellv[index];
}

////////////////////////////////////////////////////////////////////

const Array& WavelengthGrid::lambdav() const
{
    return _lambdav;
}

////////////////////////////////////////////////////////////////////

const Array& WavelengthGrid::dlambdav() const
{
    return _dlambdav;
}

////////////////////////////////////////////////////////////////////
