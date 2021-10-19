/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DEFAULTMETALLICITYCUTSPROBE_HPP
#define DEFAULTMETALLICITYCUTSPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** For each medium in the simulation that has a metallicity state variable (as requested by the
    associated material mix), DefaultMetallicityCutsProbe outputs FITS files with cuts through the
    metallicity along the coordinate planes. The field of view of each cut covers the extent of the
    spatial grid in the simulation in the relevant directions. Each cut has 1024 x 1024 pixels.

    The files are named <tt>prefix_probe_Z_XX_N.fits</tt> where XX indicates the orientation of the
    cut and N is replaced with the zero-based index of the medium in the configuration (i.e. in the
    ski file). Each file contains a single image frame.

    Note that dust components do not store the imported metallicity; for those components,
    metallicity is simply used as a multiplier to calculate the mass density. */
class DefaultMetallicityCutsProbe : public Probe
{
    /** The enumeration type indicating when probing occurs. */
    ENUM_DEF(ProbeAfter, Setup, Run)
        ENUM_VAL(ProbeAfter, Setup, "after setup")
        ENUM_VAL(ProbeAfter, Run, "after the complete simulation run")
    ENUM_END()

    ITEM_CONCRETE(DefaultMetallicityCutsProbe, Probe, "cuts of the metallicity along the coordinate planes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(DefaultMetallicityCutsProbe, "Level2&Gas&SpatialGrid")

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
