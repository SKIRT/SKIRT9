/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTTEMPERATUREPERCELLPROBE_HPP
#define DUSTTEMPERATUREPERCELLPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** DustTemperaturePerCellProbe outputs a column text file (named <tt>prefix_probe_T.dat</tt>)
    listing an indicative dust temperature for each cell in the spatial grid of the simulation.
    Specifically, the output file contains a line for each cell in the spatial grid of the
    simulation. The first column specifies the cell index, and the second column lists the
    indicative dust temperature.

    The indicative dust temperature for a spatial cell is obtained as follows. For each material
    mix of type dust present in the cell, or if applicable, for each dust population in these
    mixes, the probe calculates the equilibrium temperature that would be reached when the dust is
    embedded in the radiation field tracked by the simulation for the cell. This is achieved by
    solving the energy balance equation under LTE (local thermal equilibrium) assumptions. The
    resulting temperatures are finally averaged over the dust populations in each mix (weighed by
    the relative mass in the mix) over and all dust components present in the spatial cell (weighed
    by relative mass in the cell).

    Note that the indicative dust temperature does not correspond to a physical temperature. The
    LTE assumption is almost certainly unjustified for a relevant portion of the dust grains
    (depending on the embedding radiation field), and even when ignoring this problem, averaging
    temperatures over dust populations and dust mixes has no clear-cut physical interpretation. */
class DustTemperaturePerCellProbe : public Probe
{
    ITEM_CONCRETE(DustTemperaturePerCellProbe, Probe, "the indicative dust temperature for each spatial cell")
        ATTRIBUTE_TYPE_DISPLAYED_IF(DustTemperaturePerCellProbe, "Level2&Dust&SpatialGrid&RadiationField")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after all photon packets have been emitted and detected. */
    void probeRun() override;
};

////////////////////////////////////////////////////////////////////

#endif
