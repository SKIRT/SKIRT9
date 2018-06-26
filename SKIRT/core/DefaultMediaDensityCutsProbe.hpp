/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DEFAULTmediaDENSITYCUTSPROBE_HPP
#define DEFAULTmediaDENSITYCUTSPROBE_HPP

#include "Probe.hpp"
#include "Geometry.hpp"

////////////////////////////////////////////////////////////////////

/** DefaultmediaDensityCutsProbe outputs FITS files with cuts through the theoretical media mass
    density and the grid-discretized media mass density along the coordinate planes. Each of these
    maps contains 1024 x 1024 pixels, and covers as a field of view the total extension of the
    spatial grid in the simulation.

    The number of data files written depends on the dimension of the media system's geometry: for
    spherical symmetry only the intersection with the xy plane is written, for axial symmetry the
    intersections with the xy and xz planes are written, and for general geometries all three
    intersections are written. The difference between the theoretical mass density maps (named
    <tt>prefix_probe_trhoXX.fits</tt>) and the grid-discretized mass density maps (named
    <tt>prefix_probe_grhoXX.fits</tt>) is the following: the theoretical density is the density
    represented by the media system, i.e.\ the true dust density that would correspond to an
    infinitely fine spatial grid. The grid-discretized density maps on the other hand give the
    value of the density as read from the finite-resolution spatial grid in the simulation. A
    comparison of both sets of maps can reveal whether the chosen dust grid is suitable (in the
    ideal case, there would be no difference between both sets of maps).

    IMPLEMENTATION NOTE: the current implementation only outputs the theoretical density for a
    single geometry. The user must provide the geometry and the coordinates of the 3D
    field-of-view. Future versions will of course derive this info automatically. */
class DefaultMediaDensityCutsProbe : public Probe
{
    ITEM_CONCRETE(DefaultMediaDensityCutsProbe, Probe, "cuts of the media mass density along the coordinate axes")

    PROPERTY_ITEM(geometry, Geometry, "the geometry to be probed")
        ATTRIBUTE_DEFAULT_VALUE(geometry, "PlummerGeometry")

    PROPERTY_DOUBLE(minX, "the start point of the FOV in the X direction")
        ATTRIBUTE_QUANTITY(minX, "length")

    PROPERTY_DOUBLE(maxX, "the end point of the FOV in the X direction")
        ATTRIBUTE_QUANTITY(maxX, "length")

    PROPERTY_DOUBLE(minY, "the start point of the FOV in the Y direction")
        ATTRIBUTE_QUANTITY(minY, "length")

    PROPERTY_DOUBLE(maxY, "the end point of the FOV in the Y direction")
        ATTRIBUTE_QUANTITY(maxY, "length")

    PROPERTY_DOUBLE(minZ, "the start point of the FOV in the Z direction")
        ATTRIBUTE_QUANTITY(minZ, "length")

    PROPERTY_DOUBLE(maxZ, "the end point of the FOV in the Z direction")
        ATTRIBUTE_QUANTITY(maxZ, "length")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after setup. */
    void probeSetup() override;
};

////////////////////////////////////////////////////////////////////

#endif
