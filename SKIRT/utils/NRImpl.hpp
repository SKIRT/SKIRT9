/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NRIMPL_HPP
#define NRIMPL_HPP

#include "Array.hpp"

////////////////////////////////////////////////////////////////////

/** This namespace contains private utilities for the NR namespace. These non-template functions
    are seperated out into a seperate compilation unit (1) to avoid code-duplication and (2) to
    avoid propagating the dependencies (includes) of the implementation to all NR users. */
namespace NR_Impl
{
    /** This function calculates the normalized cdf \em Pv for the given axis grid \em xv and
        corresponding unnormalized pdf \em pv. It also normalizes the incoming pdf \em pv and
        returns the normalization factor.

        If the \em loglog flag is false, piece-wise linear behavior of both the pdf and cdf is
        assumed and regular trapezium-rule integration is used. If the \em loglog flag is true, it
        is assumed that the pdf is linear in log-log space between any two grid points (equivalent
        to power-law behavior), and the integration is performed accordingly, as described below.

        Consider the pdf values \f$p_i\f$ and \f$p_{i+1}\f$ at two consecutive grid points
        \f$x_i\f$ and \f$x_{i+1}\f$. Assuming power-law behavior, the pdf between these two grid
        points can be written as \f[ p(x) = p_i \left(\frac{x}{x_i}\right)^{\alpha_i},
        \quad\mathrm{with}\; \alpha_i = \frac{\ln p_{i+1}/\ln p_i}{\ln x_{i+1}/\ln x_i} \f]

        The area under the curve is then \f[ \int_{x_i}^{x_{i+1}} p(x)\,\mathrm{d}x =
        \int_{x_i}^{x_{i+1}} p_i \left(\frac{x}{x_i}\right)^{\alpha_i}\mathrm{d}x = p_i x_i
        \;\mathrm{gln}\left(-\alpha_i, \frac{x_{i+1}}{x_i}\right) \f] where \f$\mathrm{gln}(a,x)\f$
        is the generalized logarithm defined in the description of the SpecialFunctions::gln()
        function. */
    double cdf2(bool loglog, const Array& xv, Array& pv, Array& Pv);
}

////////////////////////////////////////////////////////////////////

#endif
