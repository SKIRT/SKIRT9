/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BOLOMETRICLUMINOSITYNORMALIZATION_HPP
#define BOLOMETRICLUMINOSITYNORMALIZATION_HPP

#include "LuminosityNormalization.hpp"

////////////////////////////////////////////////////////////////////

/** A BolometricLuminosityNormalization instance sets the normalization of a primary source by
    defining the total bolometric luminosity. */
class BolometricLuminosityNormalization : public LuminosityNormalization
{
    ITEM_CONCRETE(BolometricLuminosityNormalization, LuminosityNormalization,
                  "source normalization through the bolometric luminosity")

    PROPERTY_DOUBLE(luminosity, "the bolometric luminosity for the source")
        ATTRIBUTE_QUANTITY(luminosity, "bolluminosity")
        ATTRIBUTE_MIN_VALUE(luminosity, "]0")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the bolometric luminosity of a source with the specified %SED.
        For the present type of normalization, this function is trivial as the bolometric
        luminosity is directly configured by the user. */
    double luminosity(SED* sed) const override;
};

////////////////////////////////////////////////////////////////////

#endif
