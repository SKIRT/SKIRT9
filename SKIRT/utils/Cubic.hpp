/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CUBIC_HPP
#define CUBIC_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** This static class offers some basic functions to raise values to the third power. */
class Cubic
{
public:
    /** This function returns \f$x^3\f$. */
    static double pow3(double x) { return x * x * x; }

    /** This function returns \f$x_1^3 - x_0^3 = (x_1-x_0)(x_1^2 + x_1 x_0 + x_0^2)\f$. The second
        form is used because it is more numerically stable. */
    static double pow3(double x0, double x1) { return (x1 - x0) * (x1 * x1 + x1 * x0 + x0 * x0); }
};

////////////////////////////////////////////////////////////////////

#endif
