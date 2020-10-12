/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DEFAULTCUSTOMSTATECUTSPROBE_HPP
#define DEFAULTCUSTOMSTATECUTSPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** For each medium in the simulation that has one or more custom medium state variables (as
    requested by the associated material mix), DefaultCustomStateCutsProbe outputs a FITS file with
    cuts through the values of those custom state variables along the coordinate planes. The field
    of view of each cut covers the extent of the spatial grid in the simulation in the relevant
    directions. Each cut has 1024 x 1024 pixels.

    The files are named <tt>prefix_probe_customstate_XX_N.fits</tt> where XX indicates the
    orientation of the cut and N is replaced with the zero-based index of the medium in the
    configuration (i.e. in the ski file). Each file contains an image frame for each of the custom
    state variables of the corresponding medium, in the order in which those variables are
    requested by the associated material mix.

    \note The current implementation assumes that all custom state variables for a given medium
    component are the same physical quantity type and thus also have the same output units. In
    principle this restriction could be lifted but in that case it is unclear where to put the unit
    information in the FITS header. */
class DefaultCustomStateCutsProbe : public Probe
{
    /** The enumeration type indicating when probing occurs. */
    ENUM_DEF(ProbeAfter, Setup, Run)
        ENUM_VAL(ProbeAfter, Setup, "after setup")
        ENUM_VAL(ProbeAfter, Run, "after the complete simulation run")
    ENUM_END()

    ITEM_CONCRETE(DefaultCustomStateCutsProbe, Probe, "cuts of the custom medium state along the coordinate planes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(DefaultCustomStateCutsProbe, "Level2&CustomMediumState&SpatialGrid")

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
