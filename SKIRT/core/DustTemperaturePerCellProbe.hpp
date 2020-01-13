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

    The indicative dust temperature for a spatial cell is obtained by averaging the LTE equilibrium
    temperatures for the various dust mixes present in the cell. Note that the indicative dust
    temperature does not really correspond to a physical temperature. For more information about
    the indicative dust temperature, refer to the MediumSystem::indicativeDustTemperature()
    function. */
class DustTemperaturePerCellProbe : public Probe
{
    ITEM_CONCRETE(DustTemperaturePerCellProbe, Probe, "the indicative dust temperature for each spatial cell")
        ATTRIBUTE_TYPE_DISPLAYED_IF(DustTemperaturePerCellProbe, "Level2&Dust&SpatialGrid&RadiationField&Panchromatic")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after all photon packets have been emitted and detected. */
    void probeRun() override;
};

////////////////////////////////////////////////////////////////////

#endif
