/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PLANCKFUNCTION_HPP
#define PLANCKFUNCTION_HPP

#include "Basics.hpp"

//////////////////////////////////////////////////////////////////////

/** This simple class represents the Planck function \f[ B_\lambda(T) = \frac{2hc^2}{\lambda^5}\,
    \frac{1}{\exp\left(\dfrac{hc}{\lambda kT}\right)-1}. \f] The class has the temperature \f$T\f$
    as its only data member. */
class PlanckFunction
{
public:
    /** Basic constructor for the PlanckFunction class; just reads in and stores the temperature
        \f$T\f$. */
    explicit PlanckFunction(double T);

    /** The function returns the Planck function \f$B_\lambda(T)\f$ for an arbitrary wavelength
        \f$\lambda\f$. */
    double operator() (const double lambda) const;

    /** The function returns the temperature derivative \f$B_\lambda'(T)\f$ of the Planck function
        for an arbitrary wavelength \f$\lambda\f$. */
    double derivative(const double lambda) const;

private:
    double _T{0.};
};

//////////////////////////////////////////////////////////////////////

#endif
