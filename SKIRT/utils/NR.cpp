/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "NR.hpp"
#include "SpecialFunctions.hpp"

//////////////////////////////////////////////////////////////////////

double NR::cdf(bool loglog, Array& xv, Array& pv, Array& Pv, const Array& inxv, const Array& inpv, Range xrange)
{
    // copy the relevant portion of the axis grid
    size_t minRight = std::upper_bound(begin(inxv), end(inxv), xrange.min()) - begin(inxv);
    size_t maxRight = std::lower_bound(begin(inxv), end(inxv), xrange.max()) - begin(inxv);
    size_t n = 1 + maxRight - minRight;  // n = number of bins
    xv.resize(n+1);                      // n+1 = number of border points
    size_t i = 0;                        // i = index in target array
    xv[i++] = xrange.min();              // j = index in input array
    for (size_t j = minRight; j < maxRight; ) xv[i++] = inxv[j++];
    xv[i++] = xrange.max();

    // interpolate or copy the corresponding probability density values
    pv.resize(n+1);
    pv[0] = minRight == 0 ? 0. :
            ( loglog ? interpolateLogLog(xv[0], inxv[minRight-1], inxv[minRight], inpv[minRight-1], inpv[minRight])
                     : interpolateLinLin(xv[0], inxv[minRight-1], inxv[minRight], inpv[minRight-1], inpv[minRight]) );
    for (size_t i = 1; i < n; ++i) pv[i] = inpv[minRight+i-1];
    pv[n] = maxRight == inxv.size() ? 0. :
            ( loglog ? interpolateLogLog(xv[n], inxv[maxRight-1], inxv[maxRight], inpv[maxRight-1], inpv[maxRight])
                     : interpolateLinLin(xv[n], inxv[maxRight-1], inxv[maxRight], inpv[maxRight-1], inpv[maxRight]) );

    // calculate cumulative values corresponding to each x grid point (and any extra axis values)
    Pv.resize(n+1);                     // also sets Pv[0] to zero
    for (size_t i = 0; i!=n; ++i)
    {
        double area = 0.;
        if (!loglog || pv[i]==0 || pv[i+1]==0)
        {
            area = 0.5*(pv[i]+pv[i+1])*(xv[i+1]-xv[i]);
        }
        else
        {
            double alpha = log(pv[i+1]/pv[i]) / log(xv[i+1]/xv[i]);
            area = pv[i]*xv[i] * SpecialFunctions::gln(-alpha, xv[i+1]/xv[i]);
        }
        Pv[i+1] = Pv[i] + area;
    }

    // normalize both cumulative and regular distribution
    double norm = Pv[n];
    pv /= norm;
    Pv /= norm;
    return norm;
}


//////////////////////////////////////////////////////////////////////
