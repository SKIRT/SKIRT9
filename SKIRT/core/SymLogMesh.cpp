/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SymLogMesh.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

Array SymLogMesh::mesh() const
{
    Array tv;
    int n = numBins();

    // small number of bins
    if (n <= 2)
    {
        NR::buildLinearGrid(tv, 0., 1., n);
    }

    // larger number of bins
    else
    {
        // build the rightmost grid, without the anchor
        Array tmpv;
        int n2 = (n - 1) / 2;
        NR::buildLogGrid(tmpv, _centralBinFraction, 1., n2);

        // make room for all grid borders and copy them
        tv.resize(n + 1);
        int k = 0;
        tv[k++] = 0.;
        for (int i = n2 - 1; i >= 0; --i) tv[k++] = 0.5 - 0.5 * tmpv[i];
        if (n % 2 == 0) tv[k++] = 0.5;
        for (int i = 0; i <= n2 - 1; ++i) tv[k++] = 0.5 + 0.5 * tmpv[i];
        tv[k++] = 1.;
    }
    return tv;
}

//////////////////////////////////////////////////////////////////////
