/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PLANARDUSTTEMPERATURECUTSPROBE_HPP
#define PLANARDUSTTEMPERATURECUTSPROBE_HPP

#include "AbstractPlanarCutsProbe.hpp"

////////////////////////////////////////////////////////////////////

/** PlanarDustTemperatureCutsProbe outputs FITS files (named <tt>prefix_probe_dust_T_XX.fits</tt>)
    with cuts through an indicative dust temperature along three planes parallel to the coordinate
    planes. The offset of each cut plane from the corresponding coordinate plane can be configured
    by the user (and is zero by default). The field of view of each cut covers the extent of the
    spatial grid in the simulation in the relevant directions. The number of pixels in each
    direction can be configured by the user as well.

    The indicative dust temperature assigned to each output pixel is obtained as follows. First the
    probe determines the cell in the simulation's spatial grid "containing" the pixel according to
    its position in the cut being considered. The probe then calculates the indicative dust
    temperature for that cell by averaging the LTE equilibrium temperatures for the various dust
    mixes present in the cell. Note that the indicative dust temperature does not really correspond
    to a physical temperature. For more information about the indicative dust temperature, refer to
    the MediumSystem::indicativeDustTemperature() function. */
class PlanarDustTemperatureCutsProbe : public AbstractPlanarCutsProbe
{
    ITEM_CONCRETE(PlanarDustTemperatureCutsProbe, AbstractPlanarCutsProbe,
                  "cuts of the indicative dust temperature along planes parallel to the coordinate planes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(PlanarDustTemperatureCutsProbe,
                                    "Level2&Dust&SpatialGrid&RadiationField&Panchromatic")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after all photon packets have been emitted and detected. */
    void probeRun() override;

    /** This function outputs a FITS file with an indicative dust temperature cut in a plane
        parallel to the coordinate plane indicated by the boolean "direction" arguments \em xd, \em
        yd, and \em zd, exactly two of which must be true. The arguments \em xc, \em yc, and \em zc
        specify the position of the cuts, and the arguments \em Nx, \em Ny, and \em Nz specify the
        number of pixels in each direction. */
    static void writeDustTemperatureCut(Probe* probe, bool xd, bool yd, bool zd,
                                        double xc, double yc, double zc, int Nx, int Ny, int Nz);
};

////////////////////////////////////////////////////////////////////

#endif
