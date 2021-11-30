/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ELECTRONTEMPERATUREPERCELLPROBE_HPP
#define ELECTRONTEMPERATUREPERCELLPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** ElectronTemperaturePerCellProbe outputs a column text file (named <tt>prefix_probe_T.dat</tt>)
    listing the electron temperature defined by the input model for each cell in the spatial grid
    of the simulation. Specifically, the output file contains a line for each cell in the spatial
    grid of the simulation. The first column specifies the cell index, and the second column lists
    the electron temperature. */
class ElectronTemperaturePerCellProbe : public Probe
{
    ITEM_CONCRETE(ElectronTemperaturePerCellProbe, Probe, "the electron temperature for each spatial cell")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ElectronTemperaturePerCellProbe, "Level2&ElectronMix")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing at the end of the setup phase. */
    void probeSetup() override;
};

////////////////////////////////////////////////////////////////////

#endif
