/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LASERANGULARDISTRIBUTION_HPP
#define LASERANGULARDISTRIBUTION_HPP

#include "AngularDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** The LaserAngularDistribution class describes the angular distribition for a point source that
    emits all its radiation towards the positive Z-axis. In other words, the probability
    distribution is a Dirac delta function, \f$4\pi\delta({\bf{k}}-{\bf{e}}_z)\f$.

    The emission pattern is axisymmetric, so this geometry has a dimension of 2. */
class LaserAngularDistribution : public AngularDistribution
{
    ITEM_CONCRETE(LaserAngularDistribution, AngularDistribution, "a laser emission profile")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the angular distribution, which is 2 in this case.
        */
    int dimension() const override;

    /** This function returns the normalized probability for a given direction \f${\bf{k}} =
        (\theta,\phi)\f$. In this case, the function returns infinity if \f${\bf{k}} =
        {\bf{e}}_z\f$, or equivalently if \f$\theta=0\f$, and zero in all other cases. */
    double probabilityForDirection(Direction bfk) const override;

    /** This function generates a random direction drawn from the angular probability distribution.
        In this case, this function returns \f${\bf{k}} = {\bf{e}}_z\f$. */
    Direction generateDirection() const override;
};

////////////////////////////////////////////////////////////////////

#endif
