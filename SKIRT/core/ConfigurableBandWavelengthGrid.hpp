/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CONFIGURABLEBANDWAVELENGTHGRID_HPP
#define CONFIGURABLEBANDWAVELENGTHGRID_HPP

#include "BandWavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the ConfigurableBandWavelengthGrid class represents a wavelength grid where each
    bin is defined by a Band object. These bands can be individually configured by the user.
    Examples include the standard Johnson filters and the broadbands for actual observatories such
    as GALEX, SDSS or Herschel. See the Band class for more information.

    The Band objects are automatically sorted in order of increasing pivot wavelength, even if the
    configuration file specifies a different order. The intervals in which the transmission is
    nonzero is allowed to overlap, but no two bands in the list should have the same pivot
    wavelength (this restriction essentially disallows specifying the same band twice). */
class ConfigurableBandWavelengthGrid : public BandWavelengthGrid
{
    ITEM_CONCRETE(ConfigurableBandWavelengthGrid, BandWavelengthGrid,
                  "a wavelength grid including a configurable list of (broad)bands")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ConfigurableBandWavelengthGrid, "Level2")

        PROPERTY_ITEM_LIST(bands, Band, "the (broad)bands defining this wavelength grid")
        ATTRIBUTE_DEFAULT_VALUE(bands, "BroadBand")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function returns the list of user-configured Band objects. */
    vector<Band*> bandList() override;
};

//////////////////////////////////////////////////////////////////////

#endif
