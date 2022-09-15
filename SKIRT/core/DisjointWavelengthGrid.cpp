/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DisjointWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

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
    _borderv.resize(n + 1);
    if (n == 1)
    {
        // just a single wavelength -> form narrow bin
        _lambdaleftv[0] = _borderv[0] = _lambdav[0] * 0.999;
        _lambdarightv[0] = _borderv[1] = _lambdav[0] * 1.001;
    }
    else
    {
        if (logScale)
        {
            _lambdaleftv[0] = _borderv[0] = sqrt(_lambdav[0] * _lambdav[0] * _lambdav[0] / _lambdav[1]);
            for (size_t ell = 1; ell != n; ++ell)
            {
                _lambdarightv[ell - 1] = _lambdaleftv[ell] = _borderv[ell] = sqrt(_lambdav[ell - 1] * _lambdav[ell]);
            }
            _lambdarightv[n - 1] = _borderv[n] =
                sqrt(_lambdav[n - 1] * _lambdav[n - 1] * _lambdav[n - 1] / _lambdav[n - 2]);
        }
        else
        {
            _lambdaleftv[0] = _borderv[0] = (3. * _lambdav[0] - _lambdav[1]) / 2.;
            for (size_t ell = 1; ell != n; ++ell)
            {
                _lambdarightv[ell - 1] = _lambdaleftv[ell] = _borderv[ell] = (_lambdav[ell - 1] + _lambdav[ell]) / 2.;
            }
            _lambdarightv[n - 1] = _borderv[n] = (3. * _lambdav[n - 1] - _lambdav[n - 2]) / 2.;
        }
    }

    // verify that all bin borders are positive
    if (_lambdaleftv[0] <= 0.0) throw FATALERROR("All wavelength bin borders should be positive");

    // calculate the bin widths
    _dlambdav = _lambdarightv - _lambdaleftv;

    // setup the mapping from border bin indices to actual wavelength bin indices (see bin() function)
    _ellv.resize(n + 2);
    _ellv[0] = -1;
    for (size_t ell = 0; ell != n; ++ell) _ellv[ell + 1] = ell;
    _ellv[n + 1] = -1;
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
    _borderv.resize(2 * n);
    if (!constantWidth)
    {
        for (size_t ell = 0; ell != n; ++ell)
        {
            _borderv[2 * ell] = _lambdaleftv[ell] = _lambdav[ell] * (1. - relativeHalfWidth);
            _borderv[2 * ell + 1] = _lambdarightv[ell] = _lambdav[ell] * (1. + relativeHalfWidth);
        }
    }
    else
    {
        double delta = _lambdav[0] * relativeHalfWidth;
        for (size_t ell = 0; ell != n; ++ell)
        {
            _borderv[2 * ell] = _lambdaleftv[ell] = _lambdav[ell] - delta;
            _borderv[2 * ell + 1] = _lambdarightv[ell] = _lambdav[ell] + delta;
        }
    }

    // verify that the bins do not overlap
    if (!std::is_sorted(begin(_borderv), end(_borderv)))
        throw FATALERROR("Non-adjacent wavelength bins should not overlap");

    // calculate the bin widths
    _dlambdav = _lambdarightv - _lambdaleftv;

    // setup the mapping from border bin indices to actual wavelength bin indices (see bin() function)
    _ellv.resize(2 * n + 1);
    _ellv[0] = -1;
    for (size_t ell = 0; ell != n; ++ell)
    {
        _ellv[2 * ell + 1] = ell;
        _ellv[2 * ell + 2] = -1;  // regions between the bins are considered out of range
    }
}

////////////////////////////////////////////////////////////////////

void DisjointWavelengthGrid::setWavelengthBorders(const Array& borderv, bool logScale)
{
    // copy and sort the specified bin borders
    _borderv = borderv;
    std::sort(begin(_borderv), end(_borderv));

    // verify that there are at least two bin borders and that the smallest one is positive
    if (_borderv.size() < 2) throw FATALERROR("There must be at least two wavelength bin borders in the grid");
    if (_borderv[0] <= 0.0) throw FATALERROR("All wavelength bin borders must be positive");

    // verify that there are no duplicates
    if (std::unique(begin(_borderv), end(_borderv)) != end(_borderv))
        throw FATALERROR("There should be no duplicate wavelength bin borders in the grid");

    // make n refer to the number of bins, not the number of borders
    size_t n = _borderv.size() - 1;

    // copy the left and right bin borders
    _lambdaleftv.resize(n);
    _lambdarightv.resize(n);
    for (size_t ell = 0; ell != n; ++ell)
    {
        _lambdaleftv[ell] = _borderv[ell];
        _lambdarightv[ell] = _borderv[ell + 1];
    }

    // calculate the characteristic wavelengths
    _lambdav.resize(n);
    if (logScale)
    {
        for (size_t ell = 0; ell != n; ++ell)
        {
            _lambdav[ell] = sqrt(_lambdaleftv[ell] * _lambdarightv[ell]);
        }
    }
    else
    {
        for (size_t ell = 0; ell != n; ++ell)
        {
            _lambdav[ell] = (_lambdaleftv[ell] + _lambdarightv[ell]) / 2.;
        }
    }

    // calculate the bin widths
    _dlambdav = _lambdarightv - _lambdaleftv;

    // setup the mapping from border bin indices to actual wavelength bin indices (see bin() function)
    _ellv.resize(n + 2);
    _ellv[0] = -1;
    for (size_t ell = 0; ell != n; ++ell) _ellv[ell + 1] = ell;
    _ellv[n + 1] = -1;
}

////////////////////////////////////////////////////////////////////

void DisjointWavelengthGrid::setWavelengthSegments(const Array& bordcharv)
{
    // verify that the number of values is uneven and least three, and that the first and last values are positive
    size_t n = bordcharv.size();
    if (n < 3) throw FATALERROR("There must be at least three wavelength values in the list");
    if (n % 2 == 0) throw FATALERROR("There number of wavelength values in the list must be uneven");
    if (bordcharv[0] <= 0. || bordcharv[n - 1] <= 0.) throw FATALERROR("All wavelength bin borders must be positive");

    // copy the values into vectors that can be passed to the function that will do the actual work
    vector<double> borderv;
    vector<double> characv;
    size_t i = 0;
    while (true)
    {
        borderv.push_back(bordcharv[i++]);
        if (i == n) break;
        characv.push_back(bordcharv[i++]);
    }

    // reverse the lists if required
    if (bordcharv[0] > bordcharv[n - 1])
    {
        std::reverse(borderv.begin(), borderv.end());
        std::reverse(characv.begin(), characv.end());
    }

    // add the "characteristic wavelength" for the final segment outside the grid
    characv.push_back(0.);

    // verify the ordering
    n = borderv.size() - 1;
    for (size_t i = 0; i != n; ++i)
    {
        if (borderv[i + 1] <= borderv[i])
            throw FATALERROR("Wavelength bin borders must be in strictly monotonous order");
        if (std::isinf(characv[i])) characv[i] = 0.;  // handle zero frequencies
        if (characv[i] != 0. && (characv[i] <= borderv[i] || characv[i] >= borderv[i + 1]))
            throw FATALERROR("Characteristic wavelength must be within bin borders");
    }

    // call the function that will do the actual work
    setWavelengthSegments(borderv, characv);
}

////////////////////////////////////////////////////////////////////

void DisjointWavelengthGrid::setWavelengthSegments(const vector<double>& borderv, const vector<double>& characv)
{
    // determine the number of borders and the number of actual bins
    int numBorders = borderv.size();
    int numBins = std::count_if(characv.begin(), characv.end(), [](double c) { return c > 0.; });

    // allocate memory for all arrays
    _lambdav.resize(numBins);       // index ell
    _dlambdav.resize(numBins);      // index ell
    _lambdaleftv.resize(numBins);   // index ell
    _lambdarightv.resize(numBins);  // index ell
    _borderv.resize(numBorders);    // index k
    _ellv.resize(numBorders + 1);   // index k+1

    // copy the data into our arrays, and construct the bin index mapping
    _ellv[0] = -1;
    int ell = 0;  // ell is bin index; k is border index
    for (int k = 0; k != numBorders; ++k)
    {
        _borderv[k] = borderv[k];
        if (characv[k] > 0.)
        {
            _lambdav[ell] = characv[k];
            _dlambdav[ell] = borderv[k + 1] - borderv[k];
            _lambdaleftv[ell] = borderv[k];
            _lambdarightv[ell] = borderv[k + 1];
            _ellv[k + 1] = ell++;
        }
        else
        {
            _ellv[k + 1] = -1;
        }
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

Array DisjointWavelengthGrid::extlambdav() const
{
    int n = _lambdav.size();
    Array extv(n + 2);

    extv[0] = _lambdaleftv[0];
    for (int ell = 0; ell != n; ++ell) extv[ell + 1] = _lambdav[ell];
    extv[n + 1] = _lambdarightv[n - 1];

    return extv;
}

////////////////////////////////////////////////////////////////////

Array DisjointWavelengthGrid::extdlambdav() const
{
    int n = _dlambdav.size();
    Array extv(n + 2);

    extv[0] = 0.;
    for (int ell = 0; ell != n; ++ell) extv[ell + 1] = _dlambdav[ell];
    extv[n + 1] = 0.;

    return extv;
}

////////////////////////////////////////////////////////////////////
