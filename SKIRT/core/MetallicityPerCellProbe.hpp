/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef METALLICITYPERCELLPROBE_HPP
#define METALLICITYPERCELLPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** MetallicityPerCellProbe outputs a column text file for each medium in the simulation that has a
    metallicity state variable (as requested by the associated material mix). The files are named
    <tt>prefix_probe_Z_N.dat</tt> where N is replaced with the zero-based index of the medium in
    the configuration (i.e. in the ski file). Each file contains a line for each cell in the
    spatial grid of the simulation, and each line contains a column for the cell index and one for
    the corresponding metallicity.

    Note that dust components do not store the imported metallicity; for those components,
    metallicity is simply used as a multiplier to calculate the mass density. */
class MetallicityPerCellProbe : public Probe
{
    /** The enumeration type indicating when probing occurs. */
    ENUM_DEF(ProbeAfter, Setup, Run)
        ENUM_VAL(ProbeAfter, Setup, "after setup")
        ENUM_VAL(ProbeAfter, Run, "after the complete simulation run")
    ENUM_END()

    ITEM_CONCRETE(MetallicityPerCellProbe, Probe, "metallicity values for all spatial cells")
        ATTRIBUTE_TYPE_DISPLAYED_IF(MetallicityPerCellProbe, "Level2&Gas&SpatialGrid")

        PROPERTY_ENUM(probeAfter, ProbeAfter, "when to probe the medium state")
        ATTRIBUTE_DEFAULT_VALUE(probeAfter, "Setup")
        ATTRIBUTE_DISPLAYED_IF(probeAfter, "HasDynamicState")

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
