/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NOPOLARIZATIONPROFILE_HPP
#define NOPOLARIZATIONPROFILE_HPP

#include "PolarizationProfile.hpp"

////////////////////////////////////////////////////////////////////

/** The NoPolarizationProfile class describes the polarization profile of a source that emits
    unpolarized radiation in all directions. */
class NoPolarizationProfile : public PolarizationProfile
{
    ITEM_CONCRETE(NoPolarizationProfile, PolarizationProfile, "an unpolarized emission profile")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the polarization profile, which depends on its (lack
        of) symmetry. The unpolarized emission profile is spherically symmetric, so the function
        returns 1. */
    int dimension() const override;

    /** This function returns the Stokes vector defining the polarization state of the radiation
        emitted into the given direction \f$(\theta,\phi)\f$. For the unpolarized emission profile,
        this function returns a default-constructed StokesVector instance. */
    StokesVector polarizationForDirection(Direction bfk) const override;
};

////////////////////////////////////////////////////////////////////

#endif
