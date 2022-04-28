/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RADIATIONFIELDWAVELENGTHGRIDPROBE_HPP
#define RADIATIONFIELDWAVELENGTHGRIDPROBE_HPP

#include "SpecialtyProbe.hpp"

////////////////////////////////////////////////////////////////////

/** RadiationFieldWavelengthGridProbe outputs a column text file (named
    <tt>prefix_probe_wavelengths.dat</tt>) with details on the radiation field wavelength grid,
    i.e. the wavelength grid returned by the Configuration::radiationFieldWLG() function. For each
    wavelength bin, the file lists the characteristic wavelength, the wavelength bin width, and the
    left and right borders of the bin. */
class RadiationFieldWavelengthGridProbe : public SpecialtyProbe
{
    ITEM_CONCRETE(RadiationFieldWavelengthGridProbe, SpecialtyProbe, "wavelength grid: radiation field ")
        ATTRIBUTE_TYPE_DISPLAYED_IF(RadiationFieldWavelengthGridProbe, "Level2&Medium&SpatialGrid&RadiationField")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif
