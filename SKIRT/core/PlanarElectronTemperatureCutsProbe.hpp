/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PLANARELECTRONTEMPERATURECUTSPROBE_HPP
#define PLANARELECTRONTEMPERATURECUTSPROBE_HPP

#include "AbstractPlanarCutsProbe.hpp"

////////////////////////////////////////////////////////////////////

/** PlanarElectronTemperatureCutsProbe outputs FITS files (named
    <tt>prefix_probe_elec_T_XX.fits</tt>) with cuts through the electron temperature defined by the
    input model along three planes parallel to the coordinate planes. The offset of each cut plane
    from the corresponding coordinate plane can be configured by the user (and is zero by default).
    The field of view of each cut covers the extent of the spatial grid in the simulation in the
    relevant directions. The number of pixels in each direction can be configured by the user as
    well. */
class PlanarElectronTemperatureCutsProbe : public AbstractPlanarCutsProbe
{
    ITEM_CONCRETE(PlanarElectronTemperatureCutsProbe, AbstractPlanarCutsProbe,
                  "cuts of the electron temperature along planes parallel to the coordinate planes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(PlanarElectronTemperatureCutsProbe, "ElectronMix")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing at the end of the setup phase. */
    void probeSetup() override;

    /** This function outputs a FITS file with an electron temperature cut in a plane parallel to
        the coordinate plane indicated by the boolean "direction" arguments \em xd, \em yd, and \em
        zd, exactly two of which must be true. The arguments \em xc, \em yc, and \em zc specify the
        position of the cuts, and the arguments \em Nx, \em Ny, and \em Nz specify the number of
        pixels in each direction. */
    static void writeElectronTemperatureCut(Probe* probe, bool xd, bool yd, bool zd, double xc, double yc, double zc,
                                            int Nx, int Ny, int Nz);
};

////////////////////////////////////////////////////////////////////

#endif
