/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GASTEMPERATUREPERCELLPROBE_HPP
#define GASTEMPERATUREPERCELLPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** GasTemperaturePerCellProbe outputs a column text file (named <tt>prefix_probe_T.dat</tt>)
    listing the gas temperature defined by the input model for each cell in the spatial grid of the
    simulation. Specifically, the output file contains a line for each cell in the spatial grid of
    the simulation. The first column specifies the cell index, and the second column lists the gas
    temperature.

    In the current implementation, the probe produces output only if the simulation has been
    configured for Lyman-alpha line transfer. */
class GasTemperaturePerCellProbe : public Probe
{
    ITEM_CONCRETE(GasTemperaturePerCellProbe, Probe, "the gas temperature for each spatial cell")
        ATTRIBUTE_TYPE_DISPLAYED_IF(GasTemperaturePerCellProbe, "Lya")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing at the end of the setup phase. */
    void probeSetup() override;
};

////////////////////////////////////////////////////////////////////

#endif
