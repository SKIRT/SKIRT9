/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PLANARMEDIUMVELOCITYCUTSPROBE_HPP
#define PLANARMEDIUMVELOCITYCUTSPROBE_HPP

#include "AbstractPlanarCutsProbe.hpp"

////////////////////////////////////////////////////////////////////

/** PlanarMediumVelocityCutsProbe outputs FITS files with cuts through the spatial velocity
    distribution of the medium along three planes parallel to the coordinate planes. The offset of
    each cut plane from the corresponding coordinate plane can be configured by the user (and is
    zero by default). The field of view of each cut covers the extent of the spatial grid in the
    simulation in the relevant directions. The number of pixels in each direction can be configured
    by the user as well.

    A FITS file named <tt>prefix_probe_v_XX.fits</tt> is produced for each cut. Each FITS file
    contains three image frames respectively representing the three velocity vector components in
    the frame of the cut, i.e. the two components projected on the x and y axes of the cut and the
    component perpendicular to the cut, where positive values indicate vectors pointing away from
    the viewer. */
class PlanarMediumVelocityCutsProbe : public AbstractPlanarCutsProbe
{
    ITEM_CONCRETE(PlanarMediumVelocityCutsProbe, AbstractPlanarCutsProbe,
                  "cuts of the medium velocity along planes parallel to the coordinate planes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(PlanarMediumVelocityCutsProbe, "Level2&SpatialGrid&MediumVelocity")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after setup. */
    void probeSetup() override;

    /** This function outputs FITS files with cuts through the spatial medium velocity distribution
        along a plane parallel to the coordinate plane indicated by the boolean "direction"
        arguments \em xd, \em yd, and \em zd, exactly two of which must be true. The arguments \em
        xc, \em yc, and \em zc specify the position of the cuts, and the arguments \em Nx, \em Ny,
        and \em Nz specify the number of pixels in each direction. */
    static void writeMediumVelocityCut(Probe* probe, bool xd, bool yd, bool zd, double xc, double yc, double zc, int Nx,
                                       int Ny, int Nz);
};

////////////////////////////////////////////////////////////////////

#endif
