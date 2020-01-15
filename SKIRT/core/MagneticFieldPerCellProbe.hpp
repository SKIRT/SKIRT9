/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MAGNETICFIELDPERCELLPROBE_HPP
#define MAGNETICFIELDPERCELLPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** MagneticFieldPerCellProbe outputs a column text file (named <tt>prefix_probe_B.dat</tt>)
    listing the magnetic field for each cell in the spatial grid of the simulation.
    Specifically, the output file contains a line for each cell in the spatial grid of the
    simulation. The first column specifies the cell index, and the second, third and fourth
    column list the magnetic field components. */
class MagneticFieldPerCellProbe : public Probe
{
    ITEM_CONCRETE(MagneticFieldPerCellProbe, Probe, "the magnetic field for each spatial cell")
        ATTRIBUTE_TYPE_DISPLAYED_IF(MagneticFieldPerCellProbe, "Level3&Medium&SpatialGrid&MagneticField")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after all photon packets have been emitted and detected. */
    void probeRun() override;
};

////////////////////////////////////////////////////////////////////

#endif
