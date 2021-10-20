/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GASTEMPERATUREPERCELLPROBE_HPP
#define GASTEMPERATUREPERCELLPROBE_HPP

#include "StateProbe.hpp"

////////////////////////////////////////////////////////////////////

/** GasTemperaturePerCellProbe outputs a column text file (named <tt>prefix_probe_T.dat</tt>)
    listing the indicative gas temperature for each cell in the spatial grid of the simulation.
    Specifically, the output file contains a line for each cell in the spatial grid of the
    simulation. The first column specifies the cell index, and the second column lists the gas
    temperature. */
class GasTemperaturePerCellProbe : public StateProbe
{
    ITEM_CONCRETE(GasTemperaturePerCellProbe, StateProbe, "the gas temperature for each spatial cell")
        ATTRIBUTE_TYPE_DISPLAYED_IF(GasTemperaturePerCellProbe, "Level2&Gas&SpatialGrid")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs the probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif
