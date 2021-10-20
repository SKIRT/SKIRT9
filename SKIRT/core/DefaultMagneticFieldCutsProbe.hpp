/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DEFAULTMAGNETICFIELDCUTSPROBE_HPP
#define DEFAULTMAGNETICFIELDCUTSPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** DefaultMagneticFieldCutsProbe outputs FITS files with cuts through the spatial magnetic field
    distribution along the coordinate planes. The field of view of each cut covers the extent of
    the spatial grid in the simulation in the relevant directions. Each cut has 1024 x 1024 pixels.

    The number of data files written depends on the geometry and material contents of the media
    system. For spherical symmetry only the intersection with the xy plane is written, for axial
    symmetry the intersections with the xy and xz planes are written, and for general geometries
    all three intersections are written.

    A FITS file named <tt>prefix_probe_B_XX.fits</tt> is produced for each cut. Each FITS file
    contains three image frames respectively representing the three magnetic field vector
    components in the frame of the cut, i.e. the two components projected on the x and y axes of
    the cut and the component perpendicular to the cut, where positive values indicate vectors
    pointing away from the viewer. */
class DefaultMagneticFieldCutsProbe : public Probe
{
    ITEM_CONCRETE(DefaultMagneticFieldCutsProbe, Probe, "cuts of the magnetic field along the coordinate planes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(DefaultMagneticFieldCutsProbe, "Level3&SpatialGrid&MagneticField")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after setup. */
    void probeSetup() override;
};

////////////////////////////////////////////////////////////////////

#endif
