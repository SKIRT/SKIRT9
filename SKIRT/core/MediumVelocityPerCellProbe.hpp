/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEDIUMVELOCITYPERCELLPROBE_HPP
#define MEDIUMVELOCITYPERCELLPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** MediumVelocityPerCellProbe outputs a column text file (named <tt>prefix_probe_v.dat</tt>)
    listing the bulk velocity of the medium for each cell in the spatial grid of the simulation.
    Specifically, the output file contains a line for each cell in the spatial grid of the
    simulation. The first column specifies the cell index, and the second, third and fourth
    column list the mvelocity components. */
class MediumVelocityPerCellProbe : public Probe
{
    ITEM_CONCRETE(MediumVelocityPerCellProbe, Probe, "the medium velocity for each spatial cell")
        ATTRIBUTE_TYPE_DISPLAYED_IF(MediumVelocityPerCellProbe, "Level2&Medium&SpatialGrid&MediumVelocity")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after setup. */
    void probeSetup() override;
};

////////////////////////////////////////////////////////////////////

#endif
