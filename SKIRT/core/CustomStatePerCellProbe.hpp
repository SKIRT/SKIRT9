/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CUSTOMSTATEPERCELLPROBE_HPP
#define CUSTOMSTATEPERCELLPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** CustomStatePerCellProbe outputs a column text file for each medium in the simulation that has
    at least one custom medium state variable (as requested by the associated material mix). The
    files are named <tt>prefix_probe_customstate_N.dat</tt> where N is replaced with the zero-based
    index of the medium in the configuration (i.e. in the ski file). Each file contains a line for
    each cell in the spatial grid of the simulation, and each line contains columns representing
    the values of the custom medium state variables, in addition to the cell index. */
class CustomStatePerCellProbe : public Probe
{
    /** The enumeration type indicating when probing occurs. */
    ENUM_DEF(ProbeAfter, Setup, Run)
        ENUM_VAL(ProbeAfter, Setup, "after setup")
        ENUM_VAL(ProbeAfter, Run, "after the complete simulation run")
    ENUM_END()

    ITEM_CONCRETE(CustomStatePerCellProbe, Probe, "custom medium state variable values for all spatial cells")
        ATTRIBUTE_TYPE_DISPLAYED_IF(CustomStatePerCellProbe, "Level2&CustomMediumState&SpatialGrid")

        PROPERTY_ENUM(probeAfter, ProbeAfter, "when to probe the medium state")
        ATTRIBUTE_DEFAULT_VALUE(probeAfter, "Setup")
        ATTRIBUTE_DISPLAYED_IF(probeAfter, "DynamicState")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after setup. It produces output only if the \em
        probeAfter property is set to Setup. */
    void probeSetup() override;

    /** This function performs probing after all photon packets have been emitted and detected. It
        produces output only if the \em probeAfter property is set to Run. */
    void probeRun() override;

private:
    /** This function performs the probing; it is called from probeSetup() or probeRun() depending
        on the value of the \em probeAfter property. */
    void probe();
};

////////////////////////////////////////////////////////////////////

#endif
