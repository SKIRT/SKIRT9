/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PLANARMETALLICITYCUTSPROBE_HPP
#define PLANARMETALLICITYCUTSPROBE_HPP

#include "AbstractPlanarCutsStateProbe.hpp"

////////////////////////////////////////////////////////////////////

/** For each medium in the simulation that has a metallicity state variable (as requested by the
    associated material mix), PlanarMetallicityCutsProbe output FITS files with cuts through the
    metallicity along three planes parallel to the coordinate planes. The offset of each cut plane
    from the corresponding coordinate plane can be configured by the user (and is zero by default).
    The field of view of each cut covers the extent of the spatial grid in the simulation in the
    relevant directions. The number of pixels in each direction can be configured by the user as
    well.

    The output files are named <tt>prefix_probe_Z_XX_N.fits</tt> where XX indicates the orientation
    of the cut and N is replaced with the zero-based index of the medium in the configuration (i.e.
    in the ski file). Each file contains a single image frame.

    Note that dust components do not store the imported metallicity; for those components,
    metallicity is simply used as a multiplier to calculate the mass density. */
class PlanarMetallicityCutsProbe : public AbstractPlanarCutsStateProbe
{
    ITEM_CONCRETE(PlanarMetallicityCutsProbe, AbstractPlanarCutsStateProbe,
                  "cuts of the metallicity along planes parallel to the coordinate planes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(PlanarMetallicityCutsProbe, "Level2&Gas&SpatialGrid")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function outputs FITS files with metallicity cuts in a plane parallel to the
        coordinate plane indicated by the boolean "direction" arguments \em xd, \em yd, and \em zd,
        exactly two of which must be true. The arguments \em xc, \em yc, and \em zc specify the
        position of the cuts, and the arguments \em Nx, \em Ny, and \em Nz specify the number of
        pixels in each direction. */
    static void writeMetallicityCuts(Probe* probe, bool xd, bool yd, bool zd, double xc, double yc, double zc, int Nx,
                                     int Ny, int Nz);

protected:
    /** This function performs the probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif
