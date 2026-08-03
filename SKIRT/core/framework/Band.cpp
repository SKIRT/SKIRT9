/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Band.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

void Band::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    // obtain the data size and pointers
    _size = dataSize();
    _lambdav = wavelengthData();
    _transv = transmissionData();

    // calculate the pivot wavelength, assuming a normalized transmission curve
    double integral = 0.;
    for (size_t i = 1; i != _size; ++i)
    {
        double dlambda = _lambdav[i] - _lambdav[i - 1];
        double lambda = 0.5 * (_lambdav[i - 1] + _lambdav[i]);
        double trans = 0.5 * (_transv[i - 1] + _transv[i]);
        integral += trans * dlambda / lambda / lambda;
    }
    _pivot = 1. / sqrt(integral);

    // calculate the effective width, assuming a normalized transmission curve
    _width = 1. / *std::max_element(_transv, _transv + _size);
}

////////////////////////////////////////////////////////////////////

Range Band::wavelengthRange() const
{
    return Range(_lambdav[0], _lambdav[_size - 1]);
}

////////////////////////////////////////////////////////////////////

double Band::transmission(double wavelength) const
{
    // interpolate from the normalized transmission curve
    size_t index = std::upper_bound(_lambdav, _lambdav + _size, wavelength) - _lambdav;
    if (index == 0 || index >= _size) return 0.;
    double T =
        NR::interpolateLinLin(wavelength, _lambdav[index - 1], _lambdav[index], _transv[index - 1], _transv[index]);

    // scale to relative transmission
    return T * _width;
}

////////////////////////////////////////////////////////////////////

double Band::meanSpecificLuminosity(const Array& lambdav, const Array& pv) const
{
    // determine the range for which the SED and the transmission curve overlap (i.e. where the values may be nonzero)
    Range range(lambdav[0], lambdav[lambdav.size() - 1]);
    range.intersect(wavelengthRange());
    if (range.width() <= 0.) return 0.;

    // build a wavelength grid restricted to the overlapping range and containing all grid points from both sources
    vector<double> newlambdav;
    for (size_t i = 0; i != _size; ++i)
        if (range.contains(_lambdav[i])) newlambdav.push_back(_lambdav[i]);
    for (double lambda : lambdav)
        if (range.contains(lambda)) newlambdav.push_back(lambda);
    NR::unique(newlambdav);
    size_t newsize = newlambdav.size();
    if (newsize < 2) return 0.;

    // interpolate both sources on the new grid
    Array newpv = NR::resample<NR::interpolateLinLin>(NR::array(newlambdav), lambdav, pv);
    Array newtransv(newsize);
    for (size_t k = 0; k != newsize; ++k) newtransv[k] = transmission(newlambdav[k]);

    // perform the convolution
    double integral = 0.;
    for (size_t k = 1; k != newsize; ++k)
    {
        double dlambda = newlambdav[k] - newlambdav[k - 1];
        double Llambda = 0.5 * (newpv[k - 1] + newpv[k]);
        double trans = 0.5 * (newtransv[k - 1] + newtransv[k]);
        integral += Llambda * trans * dlambda;
    }

    // scale back to normalized transmission
    return integral / _width;
}

////////////////////////////////////////////////////////////////////

double Band::pivotWavelength() const
{
    return _pivot;
}

////////////////////////////////////////////////////////////////////

double Band::effectiveWidth() const
{
    return _width;
}

////////////////////////////////////////////////////////////////////
