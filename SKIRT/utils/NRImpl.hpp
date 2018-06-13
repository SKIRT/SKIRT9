/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NRIMPL_HPP
#define NRIMPL_HPP

#include "Array.hpp"
#include "Range.hpp"

////////////////////////////////////////////////////////////////////

/** This namespace contains private utilities for the NR namespace. These non-template functions
    are seperated out into a seperate compilation unit (1) to avoid code-duplication and (2) to
    avoid propagating the dependencies (includes) of the implementation to all NR users. */
namespace NR_Impl
{
    /** This function calculates the normalized cdf \em Pv for the given axis grid \em xv and
        corresponding unnormalized pdf \em pv. It also normalizes the incoming pdf \em pv and
        returns the normalization factor. */
    double cdf2(bool loglog, const Array& xv, Array& pv, Array& Pv);
}

////////////////////////////////////////////////////////////////////

#endif
