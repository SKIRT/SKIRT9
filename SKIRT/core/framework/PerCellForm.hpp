/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PERCELLFORM_HPP
#define PERCELLFORM_HPP

#include "SpatialGridForm.hpp"

//////////////////////////////////////////////////////////////////////

/** PerCellForm represents a spatial grid-specific probe form. Refer to the ProbeFormBridge class
    for more information about probes and forms.

    This particular form outputs a text column file listing the quantity being probed for each cell
    in the spatial grid of the simulation. Specifically, the output file contains a line for each
    cell in the spatial grid of the simulation. The first column always specifies the cell index,
    and subsequent column(s) list the quantity being probed. */
class PerCellForm : public SpatialGridForm
{
    ITEM_CONCRETE(PerCellForm, SpatialGridForm, "a text column file with values for each spatial cell")
    ITEM_END()

public:
    /** This function causes the form to output file(s) as described in the class header for the
        quantity being probed according to the information provided by the specified
        ProbeFormBridge instance. */
    void writeQuantity(const ProbeFormBridge* bridge) const override;
};

//////////////////////////////////////////////////////////////////////

#endif
