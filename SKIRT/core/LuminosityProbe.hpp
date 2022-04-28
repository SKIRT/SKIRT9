/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LUMINOSITYPROBE_HPP
#define LUMINOSITYPROBE_HPP

#include "SpecialtyWavelengthGridProbe.hpp"

////////////////////////////////////////////////////////////////////

/** LuminosityProbe outputs a text column file with the primary source luminosities on a specified
    wavelength grid (or on the default instrument wavelength grid). The file is named
    <tt>prefix_luminosities.txt</tt>. The first three columns list the characteristic wavelength of
    the wavelength bin, the total specific luminosity at that wavelength, and the total luminosity
    integrated over the bin (in both cases summed over all primary sources). Furthermore, there is
    an additional column for each primary source, listing the fractional contribution of that
    source. The columns are in the same order as the sources appear in the configuration file. If
    there is only one source, there is a single column containing one's.

    Note that, by definition, the luminosity is zero outside of the primary source wavelength range
    configured in the source system. */
class LuminosityProbe : public SpecialtyWavelengthGridProbe
{
    ITEM_CONCRETE(LuminosityProbe, SpecialtyWavelengthGridProbe, "source: luminosities of primary sources")
        ATTRIBUTE_TYPE_DISPLAYED_IF(LuminosityProbe, "Source")

    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif
