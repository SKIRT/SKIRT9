/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef METALLICITYPERCELLPROBE_HPP
#define METALLICITYPERCELLPROBE_HPP

#include "StateProbe.hpp"

////////////////////////////////////////////////////////////////////

/** MetallicityPerCellProbe outputs a column text file for each medium in the simulation that has a
    metallicity state variable (as requested by the associated material mix). The files are named
    <tt>prefix_probe_Z_N.dat</tt> where N is replaced with the zero-based index of the medium in
    the configuration (i.e. in the ski file). Each file contains a line for each cell in the
    spatial grid of the simulation, and each line contains a column for the cell index and one for
    the corresponding metallicity.

    Note that dust components do not store the imported metallicity; for those components,
    metallicity is simply used as a multiplier to calculate the mass density. */
class MetallicityPerCellProbe : public StateProbe
{
    ITEM_CONCRETE(MetallicityPerCellProbe, StateProbe, "metallicity values for all spatial cells")
        ATTRIBUTE_TYPE_DISPLAYED_IF(MetallicityPerCellProbe, "Level2&Gas&SpatialGrid")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function performs the probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif
