/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CONICALANGULARDISTRIBUTION_HPP
#define CONICALANGULARDISTRIBUTION_HPP

#include "AxAngularDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** The ConicalAngularDistribution class describes anisotropic emission that can represent a simple
    model for an AGN: the emission is isotropic within a cone around the symmetry axis (on both
    sides), and completely blocked outside this cone. The (half) opening angle \f$\Delta\f$ of the
    cone is the only free parameter. The emission pattern is axisymmetric relative to an arbitrary
    symmetry axis configured in the base class. */
class ConicalAngularDistribution : public AxAngularDistribution
{
    ITEM_CONCRETE(ConicalAngularDistribution, AxAngularDistribution, "an anisotropic conical emission profile")

        PROPERTY_DOUBLE(openingAngle, "the (half) opening angle of the cone")
        ATTRIBUTE_QUANTITY(openingAngle, "posangle")
        ATTRIBUTE_MIN_VALUE(openingAngle, "]0 deg")
        ATTRIBUTE_MAX_VALUE(openingAngle, "90 deg]")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function calculates the cosine of the opening angle for later use. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the normalized probability for a given inclination cosine
        \f$\cos\theta\f$ according to an isotropic emission within a cone with (half) opening angle
        \f$\Delta\f$. We have \f[ p({\theta}) = \frac{1}{1-\cos\Delta}\times\begin{cases}\; 1 &
        \qquad\text{if $0\le\theta<\Delta$} \\ \;0 & \qquad\text{if $\Delta\le\theta<\pi-\Delta$}
        \\ \;1 & \qquad\text{if $\pi-\Delta\le\theta\le\pi$} \end{cases} \f] */
    double probabilityForInclinationCosine(double costheta) const override;

    /** This function generates a random inclination cosine \f$\cos\theta\f$ according to an
        isotropic distribution within a cone with (half) opening angle \f$Delta\f$. This can be
        easily obtained by picking a uniform deviate \f${\cal{X}}\f$ and setting \f[ \theta =
        \begin{cases} \;\arccos\left[ 1-2{\cal{X}}(1-\cos\Delta)\right] & \qquad\text{if $0\le
        {\cal{X}} < \dfrac{\pi}{2}$} \\ \;\arccos\left[ 1-2\cos\Delta - 2{\cal{X}}(1-\cos\Delta)
        \right] & \qquad\text{if $\dfrac12 \le {\cal{X}} < 1$} \end{cases} \f] */
    double generateInclinationCosine() const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    double _cosDelta;
};

////////////////////////////////////////////////////////////////////

#endif
