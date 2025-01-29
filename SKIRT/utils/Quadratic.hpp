/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef QUADRATIC_HPP
#define QUADRATIC_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** This static class offers functions that solve quadratic equations using a numerically stable
    technique. Specifically, consider a quadratic equation of the form \f[x^2+2bx+c=0.\f] If \f$b^2
    \ge c\f$, this equation has solutions described by \f[\begin{aligned} x_1 &= -b - \sqrt{b^2-c}
    \\ x_2 &= -b + \sqrt{b^2-c} \\ s_1s_2&=c.\end{aligned}\f] To avoid loss of significance in case
    the solutions have a different order of magnitude, the functions in this class use the first
    and third equations if \f$b>0\f$ and the second and third equations otherwise. */
class Quadratic
{
public:
    /** This function determines the solutions of \f$x^2 + 2bx + c = 0\f$. If there are two
        distinct real solutions, they are stored in the arguments x1 and x2. Otherwise, i.e. if
        there are no solutions or there is just one real solution, x1 and x2 remain unchanged. */
    static void distinctSolutions(double b, double c, double& x1, double& x2)
    {
        if (b * b > c)  // if discriminant is strictly positive, there are two distinct real solutions
        {
            if (b > 0)  // x1 is always negative
            {
                x1 = -b - sqrt(b * b - c);
                x2 = c / x1;
            }
            else  // x2 is always positive
            {
                x2 = -b + sqrt(b * b - c);
                x1 = c / x2;
            }
        }
    }

    /** This function returns the smallest positive solution of \f$x^2 + 2bx + c = 0\f$, or zero
        if there is no positive solution. */
    static double smallestPositiveSolution(double b, double c)
    {
        if (b * b > c)  // if discriminant is negative, there are no real solutions
        {
            if (b > 0.)  // x1 is always negative; x2 is positive only if c<0
            {
                if (c < 0.)
                {
                    double x1 = -b - sqrt(b * b - c);
                    return c / x1;
                }
            }
            else  // x2 is always positive; x1 is positive only if c>0
            {
                double x2 = -b + sqrt(b * b - c);
                if (c > 0.)
                {
                    double x1 = c / x2;
                    if (x1 < x2) return x1;
                }
                return x2;
            }
        }
        return 0.;
    }
};

////////////////////////////////////////////////////////////////////

#endif
