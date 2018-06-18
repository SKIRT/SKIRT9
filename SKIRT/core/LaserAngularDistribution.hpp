/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LASERANGULARDISTRIBUTION_HPP
#define LASERANGULARDISTRIBUTION_HPP

#include "AxAngularDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** The LaserAngularDistribution class describes the axisymmetric angular emission distribution for
    a point source that emits all its radiation in a single direction, namely towards the positive
    symmetry axis configured in the base class. In other words, the probability distribution is a
    Dirac delta function, \f$4\pi\delta({\bf{k}}-{\bf{e}}_\mathrm{sym})\f$. */
class LaserAngularDistribution : public AxAngularDistribution
{
    ITEM_CONCRETE(LaserAngularDistribution, AxAngularDistribution, "a laser emission profile")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the normalized probability for a given inclination cosine
        \f$\cos\theta\f$ relative to the symmetry axis. In this case, the function returns infinity
        if \f$\theta=0\f$, or equivalently \f$\cos\theta=1\f$, and zero in all other cases. */
    double probabilityForInclinationCosine(double costheta) const override;

    /** This function generates a random inclination cosine relative to the symmetry axis, drawn
        from the angular probability distribution. In this case, the function returns
        \f$\cos\theta=1\f$, which is equivalent to \f$\theta=0\f$. */
    double generateInclinationCosine() const override;
};

////////////////////////////////////////////////////////////////////

#endif
