/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DEFAULTCUSTOMSTATECUTSPROBE_HPP
#define DEFAULTCUSTOMSTATECUTSPROBE_HPP

#include "StateProbe.hpp"

////////////////////////////////////////////////////////////////////

/** For each medium in the simulation that has one or more custom medium state variables (as
    requested by the associated material mix), DefaultCustomStateCutsProbe outputs a FITS file with
    cuts through the values of those custom state variables along the coordinate planes. The field
    of view of each cut covers the extent of the spatial grid in the simulation in the relevant
    directions. Each cut has 1024 x 1024 pixels.

    The files are named <tt>prefix_probe_customstate_XX_N.fits</tt> where XX indicates the
    orientation of the cut and N is replaced with the zero-based index of the medium in the
    configuration (i.e. in the ski file). Each file contains an image frame for each of the custom
    state variables of the corresponding medium, in the order in which those variables are
    requested by the associated material mix.

    \note The current implementation assumes that all custom state variables for a given medium
    component are the same physical quantity type and thus also have the same output units. In
    principle this restriction could be lifted but in that case it is unclear where to put the unit
    information in the FITS header. */
class DefaultCustomStateCutsProbe : public StateProbe
{
    ITEM_CONCRETE(DefaultCustomStateCutsProbe, StateProbe,
                  "cuts of the custom medium state along the coordinate planes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(DefaultCustomStateCutsProbe, "Level2&CustomMediumState&SpatialGrid")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function performs the probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif
