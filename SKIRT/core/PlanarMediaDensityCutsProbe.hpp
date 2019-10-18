/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PLANARMEDIADENSITYCUTSPROBE_HPP
#define PLANARMEDIADENSITYCUTSPROBE_HPP

#include "AbstractPlanarCutsProbe.hpp"

////////////////////////////////////////////////////////////////////

/** PlanarMediaDensityCutsProbe outputs FITS files with cuts through the true input media density
    and the grid-discretized media density along three planes parallel to the coordinate planes.
    The offset of each cut plane from the corresponding coordinate plane can be configured by the
    user (and is zero by default). The field of view of each cut covers the extent of the spatial
    grid in the simulation in the relevant directions. The number of pixels in each direction can
    be configured by the user as well.

    A separate set of files is produced for each material type (dust, electrons, or gas). The files
    for dust provide the dust mass density, and those for electrons and gas provide the electron or
    hydrogen number density. Multiple media components containing the same material type are
    combined, regardless of ordering.

    The difference between the true input density maps (named <tt>prefix_probe_MM_t_XX.fits</tt>)
    and the grid-discretized density maps (named <tt>prefix_probe_MM_g_XX.fits</tt>) is the
    following: the true input density is the density represented by the input model configured in
    the media system, i.e.\ the true density that would correspond to an infinitely fine spatial
    grid. The grid-discretized density maps on the other hand give the value of the density as read
    from the finite-resolution spatial grid in the simulation. A comparison of both sets of maps
    can reveal whether the configured spatial grid is suitable (in the ideal case, there would be
    no difference between both sets of maps). */
class PlanarMediaDensityCutsProbe : public AbstractPlanarCutsProbe
{
    ITEM_CONCRETE(PlanarMediaDensityCutsProbe, AbstractPlanarCutsProbe,
                  "cuts of the media densities along planes parallel to the coordinate planes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(PlanarMediaDensityCutsProbe, "Medium&SpatialGrid")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after setup. */
    void probeSetup() override;

    /** This function outputs FITS files with the theoretical and grid density for each material
        type in the coordinate planes (xy, xz, or yz) indicated by the boolean "direction"
        arguments \em xd, \em yd, and \em zd, exactly two of which must be true. The arguments \em
        xc, \em yc, and \em zc specify the position of the cuts, and the arguments \em Nx, \em Ny,
        and \em Nz specify the number of pixels in each direction. */
    static void writeMediaDensityCuts(Probe* probe, bool xd, bool yd, bool zd,
                                      double xc, double yc, double zc, int Nx, int Ny, int Nz);
};

////////////////////////////////////////////////////////////////////

#endif
