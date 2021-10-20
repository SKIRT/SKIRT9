/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PLANARGASTEMPERATURECUTSPROBE_HPP
#define PLANARGASTEMPERATURECUTSPROBE_HPP

#include "AbstractPlanarCutsStateProbe.hpp"

////////////////////////////////////////////////////////////////////

/** PlanarGasTemperatureCutsProbe outputs FITS files (named <tt>prefix_probe_dust_T_XX.fits</tt>)
    with cuts through the indicative gas temperature along three planes parallel to the coordinate
    planes. The offset of each cut plane from the corresponding coordinate plane can be configured
    by the user (and is zero by default). The field of view of each cut covers the extent of the
    spatial grid in the simulation in the relevant directions. The number of pixels in each
    direction can be configured by the user as well. */
class PlanarGasTemperatureCutsProbe : public AbstractPlanarCutsStateProbe
{
    ITEM_CONCRETE(PlanarGasTemperatureCutsProbe, AbstractPlanarCutsStateProbe,
                  "cuts of the gas temperature along planes parallel to the coordinate planes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(PlanarGasTemperatureCutsProbe, "Level2&Gas&SpatialGrid")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function outputs a FITS file with a gas temperature cut in a plane
        parallel to the coordinate plane indicated by the boolean "direction" arguments \em xd, \em
        yd, and \em zd, exactly two of which must be true. The arguments \em xc, \em yc, and \em zc
        specify the position of the cuts, and the arguments \em Nx, \em Ny, and \em Nz specify the
        number of pixels in each direction. */
    static void writeGasTemperatureCut(Probe* probe, bool xd, bool yd, bool zd, double xc, double yc, double zc, int Nx,
                                       int Ny, int Nz);

protected:
    /** This function performs the probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif
