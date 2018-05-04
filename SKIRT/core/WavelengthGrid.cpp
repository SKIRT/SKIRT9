/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "WavelengthGrid.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

void WavelengthGrid::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    // the subclass should have added wavelengths
    if (!_lambdav.size()) throw FATALERROR("Wavelength grid should have been initialized by subclass");
}

////////////////////////////////////////////////////////////////////

void WavelengthGrid::setWavelengths(const Array& lambdav)
{
    // copy and sort the specified representative wavelengths
    _lambdav = lambdav;
    std::sort(begin(_lambdav), end(_lambdav));

    // verify that there is at least one wavelength and that the smallest one is positive
    if (!_lambdav.size()) throw FATALERROR("There must be at least one wavelength in the grid");
    if (_lambdav[0] <= 0.0) throw FATALERROR("All wavelengths should be positive");

    // verify that there are no duplicates
    if (std::unique(begin(_lambdav), end(_lambdav)) != end(_lambdav))
        throw FATALERROR("There should be no duplicate wavelengths in the grid");

    // calculate the bin borders
    size_t n = _lambdav.size();
    _blambdav.resize(n+1);
    _blambdav[0] = _lambdav[0] * 0.999;
    for (size_t ell=1; ell!=n; ++ell) _blambdav[ell] = sqrt(_lambdav[ell-1]*_lambdav[ell]);
    _blambdav[n] = _lambdav[n-1] * 1.001;

    // calculate the bin widths
    for (size_t ell=0; ell!=n; ++ell) _dlambdav[ell] = _blambdav[ell+1] - _blambdav[ell];
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

double WavelengthGrid::lambdamin(int ell) const
{
    return _blambdav[ell];
}

////////////////////////////////////////////////////////////////////

double WavelengthGrid::lambdamax(int ell) const
{
    return _blambdav[ell+1];
}

////////////////////////////////////////////////////////////////////

int WavelengthGrid::ell(double lambda) const
{
    return NR::locateFail(_blambdav, lambda);
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

const Array& WavelengthGrid::blambdav() const
{
    return _blambdav;
}

////////////////////////////////////////////////////////////////////
