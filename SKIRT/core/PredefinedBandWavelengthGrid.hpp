/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PREDEFINEDBANDWAVELENGTHGRID_HPP
#define PREDEFINEDBANDWAVELENGTHGRID_HPP

#include "BandWavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the PredefinedBandWavelengthGrid class represents a wavelength grid where each
    bin is defined by a Band object. The class offers a list of predefined bands comprising the
    GALEX, SDSS, 2MASS, WISE and HERSCHEL broadbands. Each of these sets can be included or excluded
    as a whole through configuration flags. */
class PredefinedBandWavelengthGrid : public BandWavelengthGrid
{
    ITEM_CONCRETE(PredefinedBandWavelengthGrid, BandWavelengthGrid,
                  "a wavelength grid including a predefined list of (broad)bands")

        PROPERTY_BOOL(includeGALEX, "include GALEX FUV and NUV bands")
        ATTRIBUTE_DEFAULT_VALUE(includeGALEX, "true")

        PROPERTY_BOOL(includeSDSS, "include SDSS ugriz bands")
        ATTRIBUTE_DEFAULT_VALUE(includeSDSS, "true")

        PROPERTY_BOOL(include2MASS, "include 2MASS J, H and Ks bands")
        ATTRIBUTE_DEFAULT_VALUE(include2MASS, "true")

        PROPERTY_BOOL(includeWISE, "include WISE W1, W2, W3 and W4 bands")
        ATTRIBUTE_DEFAULT_VALUE(includeWISE, "DustEmission:true;false")

        PROPERTY_BOOL(includeHERSCHEL, "include HERSCHEL PACS 70,100,160 and SPIRE 250,350,500 bands")
        ATTRIBUTE_DEFAULT_VALUE(includeHERSCHEL, "DustEmission:true;false")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function returns the list of predefined Band objects, depending on the user-configured
        flags. */
    vector<Band*> bandList() override;
};

//////////////////////////////////////////////////////////////////////

#endif
