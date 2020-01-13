/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef POLARIZATIONPROFILE_HPP
#define POLARIZATIONPROFILE_HPP

#include "PolarizationProfileInterface.hpp"
#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** An instance of a PolarizationProfile subclass describes the polarization profile of the radiation
    emitted by a point source into a given direction \f$(\theta,\phi)\f$. This class inherits the
    PolarizationProfileInterface interface offering the polarizationForDirection() function which
    much be implemented by a subclass. */
class PolarizationProfile : public SimulationItem, public PolarizationProfileInterface
{
    ITEM_ABSTRACT(PolarizationProfile, SimulationItem, "a polarized emission profile")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the polarization profile, which depends on its (lack
        of) symmetry. */
    virtual int dimension() const = 0;
};

////////////////////////////////////////////////////////////////////

#endif
