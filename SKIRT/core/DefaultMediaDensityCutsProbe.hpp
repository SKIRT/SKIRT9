/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DEFAULTMEDIADENSITYCUTSPROBE_HPP
#define DEFAULTMEDIADENSITYCUTSPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** DefaultMediaDensityCutsProbe outputs FITS files with cuts through the true input media density
    and the grid-discretized media density along the coordinate planes. The field of view of each
    cut covers the extent of the spatial grid in the simulation in the relevant directions. Each
    cut has 1024 x 1024 pixels.

    The number of data files written depends on the geometry and material contents of the media
    system. For spherical symmetry only the intersection with the xy plane is written, for axial
    symmetry the intersections with the xy and xz planes are written, and for general geometries
    all three intersections are written. Also, a separate set of files is produced for each
    material type (dust, electrons, or gas). The files for dust provide the dust mass density, and
    those for electrons and gas provide the electron or hydrogen number density. Multiple media
    components containing the same material type are combined, regardless of ordering.

    The difference between the true input density maps (named <tt>prefix_probe_MM_t_XX.fits</tt>)
    and the grid-discretized density maps (named <tt>prefix_probe_MM_g_XX.fits</tt>) is the
    following: the true input density is the density represented by the input model configured in
    the media system, i.e.\ the true density that would correspond to an infinitely fine spatial
    grid. The grid-discretized density maps on the other hand give the value of the density as read
    from the finite-resolution spatial grid in the simulation. A comparison of both sets of maps
    can reveal whether the configured spatial grid is suitable (in the ideal case, there would be
    no difference between both sets of maps). */
class DefaultMediaDensityCutsProbe : public Probe
{
    /** The enumeration type indicating when probing occurs. */
    ENUM_DEF(ProbeAfter, Setup, Run)
        ENUM_VAL(ProbeAfter, Setup, "after setup")
        ENUM_VAL(ProbeAfter, Run, "after the complete simulation run")
    ENUM_END()

    ITEM_CONCRETE(DefaultMediaDensityCutsProbe, Probe, "cuts of the media densities along the coordinate planes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(DefaultMediaDensityCutsProbe, "Medium&SpatialGrid")

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
