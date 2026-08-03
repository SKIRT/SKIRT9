/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTEMISSIONWAVELENGTHGRIDPROBE_HPP
#define DUSTEMISSIONWAVELENGTHGRIDPROBE_HPP

#include "SpecialtyProbe.hpp"

////////////////////////////////////////////////////////////////////

/** DustEmissionWavelengthGridProbe outputs a column text file (named
    <tt>prefix_probe_wavelengths.dat</tt>) with details on the dust emission wavelength grid, i.e.
    the wavelength grid returned by the Configuration::dustEmissionWLG() function. For each
    wavelength bin, the file lists the characteristic wavelength, the wavelength bin width, and the
    left and right borders of the bin. */
class DustEmissionWavelengthGridProbe : public SpecialtyProbe
{
    ITEM_CONCRETE(DustEmissionWavelengthGridProbe, SpecialtyProbe, "wavelength grid: dust emission")
        ATTRIBUTE_TYPE_DISPLAYED_IF(DustEmissionWavelengthGridProbe, "Level2&DustMix&DustEmission")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif
