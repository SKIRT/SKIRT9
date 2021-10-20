/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CUSTOMSTATEPERCELLPROBE_HPP
#define CUSTOMSTATEPERCELLPROBE_HPP

#include "StateProbe.hpp"

////////////////////////////////////////////////////////////////////

/** CustomStatePerCellProbe outputs a column text file for each medium in the simulation that has
    at least one custom medium state variable (as requested by the associated material mix). The
    files are named <tt>prefix_probe_customstate_N.dat</tt> where N is replaced with the zero-based
    index of the medium in the configuration (i.e. in the ski file). Each file contains a line for
    each cell in the spatial grid of the simulation, and each line contains columns representing
    the values of the custom medium state variables, in addition to the cell index. */
class CustomStatePerCellProbe : public StateProbe
{
    ITEM_CONCRETE(CustomStatePerCellProbe, StateProbe, "custom medium state variable values for all spatial cells")
        ATTRIBUTE_TYPE_DISPLAYED_IF(CustomStatePerCellProbe, "Level2&CustomMediumState&SpatialGrid")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function performs the probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif
