/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TODDLERSLINESED_HPP
#define TODDLERSLINESED_HPP

#include "FamilySED.hpp"

////////////////////////////////////////////////////////////////////

/** ToddlersLineSED is a class representing the family of 
    Toddlers Hii region template SEDs, parameterized on age, metallicity,
    star formation efficiency, and cloud number density and scaled by particle mass.

    The line and the continuum emission are currently separated, this class is 
    used for the line emission in the SED. Most of the lines here are from the 
    Cloudy Hii regions line list.
    The intrinsic line resolution for all lines (R = $\lambda / \Delta \lambda$) 
    has been fixed to 1e5.

    TODDLERS = Time evolution of Dust Diagnostics and Line Emission from Regions
    containing young Stars.
 */
class ToddlersLineSED : public FamilySED
{
    ITEM_CONCRETE(ToddlersLineSED, FamilySED, "an Hii region line emission SED from the TODDLERS model")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ToddlersLineSED, "Level2")

        PROPERTY_DOUBLE(age, "Hii region age")
        ATTRIBUTE_QUANTITY(age, "time")
        ATTRIBUTE_MIN_VALUE(age, "[0.1 Myr")
        ATTRIBUTE_MAX_VALUE(age, "30 Myr]")
        ATTRIBUTE_DEFAULT_VALUE(age, "5.0 Myr")

        PROPERTY_DOUBLE(metallicity, "Hii region metallicity")
        ATTRIBUTE_MIN_VALUE(metallicity, "[0.001")
        ATTRIBUTE_MAX_VALUE(metallicity, "0.04]")
        ATTRIBUTE_DEFAULT_VALUE(metallicity, "0.014")

        PROPERTY_DOUBLE(SFE, "Initial star formation efficiency M*/Mtot ")
        ATTRIBUTE_MIN_VALUE(SFE, "[.025")
        ATTRIBUTE_MAX_VALUE(SFE, ".15]")
        ATTRIBUTE_DEFAULT_VALUE(SFE, ".05")

        PROPERTY_DOUBLE(cloudNumDensity, "Cloud number density")
        ATTRIBUTE_QUANTITY(cloudNumDensity, "numbervolumedensity")
        ATTRIBUTE_MIN_VALUE(cloudNumDensity, "[10")
        ATTRIBUTE_MAX_VALUE(cloudNumDensity, "1280]")
        ATTRIBUTE_DEFAULT_VALUE(cloudNumDensity, "500")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function returns a newly created SEDFamily object (which is already hooked into the
        simulation item hierachy so it will be automatically deleted) and stores the parameters for
        the specific %SED configured by the user in the specified array. */
    const SEDFamily* getFamilyAndParameters(Array& parameters) override;
};

////////////////////////////////////////////////////////////////////

#endif
