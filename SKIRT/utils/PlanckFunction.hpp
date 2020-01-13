/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PLANCKFUNCTION_HPP
#define PLANCKFUNCTION_HPP

#include "Array.hpp"
#include "Range.hpp"

//////////////////////////////////////////////////////////////////////

/** This class represents the Planck function \f[ B_\lambda(\lambda,T) = \frac{2hc^2}{\lambda^5}\,
    \frac{1}{\exp\left(\dfrac{hc}{\lambda kT}\right)-1}. \f] The temperature \f$T\f$ is specified
    in the constructor. */
class PlanckFunction
{
public:
    /** The constructor for the PlanckFunction class accepts the temperature \f$T\f$ and
        precalculates some constants for later use. */
    explicit PlanckFunction(double T);

    /** This function returns the value of the Planck function \f$B_\lambda(\lambda,T)\f$ for a
        given wavelength \f$\lambda\f$ and for the temperature \f$T\f$ specified in the
        constructor. */
    double value(double lambda) const;

    /** This function call operator returns the value of the Planck function
        \f$B_\lambda(\lambda,T)\f$ for a given wavelength \f$\lambda\f$ and for the temperature
        \f$T\f$ specified in the constructor. It is equivalent to the value() function. */
    double operator()(double lambda) const { return value(lambda); }

    /** This function constructs a tabulated normalized probability density function (pdf) and the
        corresponding normalized cumulative distribution function (cdf) for the Planck function
        with temperature \f$T\f$ (specified in the constructor) within the given wavelength range.
        As the basis for the tabulation, the function constructs a logarithmic wavelength grid with
        approximately 1000 wavelength points for each order of magitude in the wavelength range.
        The integration to determine the cdf is performed using a mechanism that assumes that the
        pdf behaves as a power-law between any two grid points. */
    double cdf(Array& lambdav, Array& pv, Array& Pv, Range lambdaRange) const;

private:
    double _T{0};
    double _f1{0};
    double _f2{0};
};

//////////////////////////////////////////////////////////////////////

#endif
