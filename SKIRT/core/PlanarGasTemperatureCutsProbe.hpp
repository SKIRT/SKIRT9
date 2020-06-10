/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PLANARGASTEMPERATURECUTSPROBE_HPP
#define PLANARGASTEMPERATURECUTSPROBE_HPP

#include "AbstractPlanarCutsProbe.hpp"

////////////////////////////////////////////////////////////////////

/** PlanarGasTemperatureCutsProbe outputs FITS files (named <tt>prefix_probe_dust_T_XX.fits</tt>)
    with cuts through the gas temperature defined by the input model along three planes parallel to
    the coordinate planes. The offset of each cut plane from the corresponding coordinate plane can
    be configured by the user (and is zero by default). The field of view of each cut covers the
    extent of the spatial grid in the simulation in the relevant directions. The number of pixels
    in each direction can be configured by the user as well.

    In the current implementation, the probe produces output only if the simulation has been
    configured for Lyman-alpha line transfer. */
class PlanarGasTemperatureCutsProbe : public AbstractPlanarCutsProbe
{
    ITEM_CONCRETE(PlanarGasTemperatureCutsProbe, AbstractPlanarCutsProbe,
                  "cuts of the gas temperature along planes parallel to the coordinate planes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(PlanarGasTemperatureCutsProbe, "Lya")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing at the end of the setup phase. */
    void probeSetup() override;

    /** This function outputs a FITS file with a gas temperature cut in a plane
        parallel to the coordinate plane indicated by the boolean "direction" arguments \em xd, \em
        yd, and \em zd, exactly two of which must be true. The arguments \em xc, \em yc, and \em zc
        specify the position of the cuts, and the arguments \em Nx, \em Ny, and \em Nz specify the
        number of pixels in each direction. */
    static void writeGasTemperatureCut(Probe* probe, bool xd, bool yd, bool zd, double xc, double yc, double zc, int Nx,
                                       int Ny, int Nz);
};

////////////////////////////////////////////////////////////////////

#endif
