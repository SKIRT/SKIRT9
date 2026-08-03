/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ISOTROPICANGULARDISTRIBUTION_HPP
#define ISOTROPICANGULARDISTRIBUTION_HPP

#include "AngularDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** The IsotropicAngularDistribution describes a uniform (isotropic) angular probability
    distribution that serves to describe the emission of a point source. */
class IsotropicAngularDistribution : public AngularDistribution
{
    ITEM_CONCRETE(IsotropicAngularDistribution, AngularDistribution, "an isotropic emission profile")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the angular distribution, which depends on its (lack
        of) symmetry. The isotropic distribution is spherically symmetric, so the function returns
        1. */
    int dimension() const override;

    /** This function returns the probability for the given direction, which for this isotropic
        distribution is always 1. */
    double probabilityForDirection(Direction bfk) const override;

    /** This function generates a random direction drawn from the probability distribution. For the
        isotropic distribution this is uniform on the unit sphere. */
    Direction generateDirection() const override;
};

////////////////////////////////////////////////////////////////////

#endif
