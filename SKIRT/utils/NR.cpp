/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "NR.hpp"
#include "SpecialFunctions.hpp"

////////////////////////////////////////////////////////////////////

void NR::reverse(Array& xv, size_t blocksize)
{
    auto first = begin(xv);
    auto last = end(xv) - blocksize;
    while (first < last)
    {
        for (size_t i = 0; i != blocksize; ++i) std::iter_swap(first + i, last + i);
        first += blocksize;
        last -= blocksize;
    }
}

////////////////////////////////////////////////////////////////////

double NR::cdf2(bool loglog, const Array& xv, Array& pv, Array& Pv)
{
    size_t n = xv.size() - 1;

    // calculate cumulative values corresponding to each x grid point (and any extra axis values)
    Pv.resize(n + 1);  // also sets Pv[0] to zero
    for (size_t i = 0; i != n; ++i)
    {
        double area = 0.;
        if (!loglog)
        {
            area = 0.5 * (pv[i] + pv[i + 1]) * (xv[i + 1] - xv[i]);
        }
        else if (pv[i] > 0 && pv[i + 1] > 0)
        {
            double alpha = log(pv[i + 1] / pv[i]) / log(xv[i + 1] / xv[i]);
            area = pv[i] * xv[i] * SpecialFunctions::gln(-alpha, xv[i + 1] / xv[i]);
        }
        Pv[i + 1] = Pv[i] + area;
    }

    // normalize both cumulative and regular distribution
    double norm = Pv[n];
    if (norm > 0.)
    {
        pv /= norm;
        Pv /= norm;
    }
    return norm;
}

////////////////////////////////////////////////////////////////////
