/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PLANARRADIATIONFIELDCUTSPROBE_HPP
#define PLANARRADIATIONFIELDCUTSPROBE_HPP

#include "AbstractPlanarCutsProbe.hpp"

////////////////////////////////////////////////////////////////////

/** PlanarRadiationFieldCutsProbe outputs FITS files (named <tt>prefix_probe_J_XX.fits</tt>) with
    cuts through the mean radiation field intensity along three planes parallel to the coordinate
    planes. The offset of each cut plane from the corresponding coordinate plane can be configured
    by the user (and is zero by default). The field of view of each cut covers the extent of the
    spatial grid in the simulation in the relevant directions. The number of pixels in each
    direction can be configured by the user as well.

    Each of the output files actually contains a datacube with a map for each bin in the wavelength
    grid returned by the Configuration::radiationFieldWLG() function. The probe offers an option to
    output a separate text column file with details on the radiation field wavelength grid. For
    each wavelength bin, the file lists the characteristic wavelength, the wavelength bin width,
    and the left and right borders of the bin. */
class PlanarRadiationFieldCutsProbe : public AbstractPlanarCutsProbe
{
    ITEM_CONCRETE(PlanarRadiationFieldCutsProbe, AbstractPlanarCutsProbe,
                  "cuts of the mean radiation field intensity along planes parallel to the coordinate planes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(PlanarRadiationFieldCutsProbe, "Level2&SpatialGrid&RadiationField")

        PROPERTY_BOOL(writeWavelengthGrid, "output a text file with the radiation field wavelength grid")
        ATTRIBUTE_DEFAULT_VALUE(writeWavelengthGrid, "false")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after all photon packets have been emitted and detected. */
    void probeRun() override;

    /** This function outputs a FITS file with a mean radiation field intensity cut in a plane
        parallel to the coordinate plane indicated by the boolean "direction" arguments \em xd, \em
        yd, and \em zd, exactly two of which must be true. The arguments \em xc, \em yc, and \em zc
        specify the position of the cuts, and the arguments \em Nx, \em Ny, and \em Nz specify the
        number of pixels in each direction. */
    static void writeRadiationFieldCut(Probe* probe, bool xd, bool yd, bool zd, double xc, double yc, double zc, int Nx,
                                       int Ny, int Nz);
};

////////////////////////////////////////////////////////////////////

#endif
