/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DEFAULTMEDIUMVELOCITYCUTSPROBE_HPP
#define DEFAULTMEDIUMVELOCITYCUTSPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** DefaultMediumVelocityCutsProbe outputs FITS files with cuts through the spatial velocity
    distribution of the medium along the coordinate planes. The field of view of each cut covers
    the extent of the spatial grid in the simulation in the relevant directions. Each cut has 1024
    x 1024 pixels.

    The number of data files written depends on the geometry and material contents of the media
    system. For spherical symmetry only the intersection with the xy plane is written, for axial
    symmetry the intersections with the xy and xz planes are written, and for general geometries
    all three intersections are written.

    A FITS file named <tt>prefix_probe_v_XX.fits</tt> is produced for each cut. Each FITS file
    contains three image frames respectively representing the three velocity vector components in
    the frame of the cut, i.e. the two components projected on the x and y axes of the cut and the
    component perpendicular to the cut, where positive values indicate vectors pointing away from
    the viewer. */
class DefaultMediumVelocityCutsProbe : public Probe
{
    ITEM_CONCRETE(DefaultMediumVelocityCutsProbe, Probe, "cuts of the medium velocity along the coordinate planes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(DefaultMediumVelocityCutsProbe, "Level2&SpatialGrid&MediumVelocity")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after setup. */
    void probeSetup() override;
};

////////////////////////////////////////////////////////////////////

#endif
