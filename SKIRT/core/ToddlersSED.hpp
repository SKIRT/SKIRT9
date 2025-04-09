/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TODDLERSSED_HPP
#define TODDLERSSED_HPP

#include "FamilySED.hpp"

////////////////////////////////////////////////////////////////////

/** A ToddlersSED class instance represents a single star-forming region %SED taken from the
    TODDLERS model with SFR-normalized spectra (integrated over the first 10 Myr) for Starburst99
    stellar populations with Kroupa IMF (0.1-100 \f$\mathrm{M}_\odot\f$) and single star evolution.
    The %SED has the highest available wavelength resolution over a range from \f$0.01\f$ to
    \f$3000~\mu\mathrm{m}\f$ and dust processing is always included. It is parameterized on
    metallicity, star formation efficiency, and cloud number density. See the ToddlersSEDFamily
    class for more information. */
class ToddlersSED : public FamilySED
{
    ITEM_CONCRETE(ToddlersSED, FamilySED, "a Toddlers SED for emission from star-forming regions")

        PROPERTY_DOUBLE(metallicity, "metallicity")
        ATTRIBUTE_MIN_VALUE(metallicity, "[0.001")
        ATTRIBUTE_MAX_VALUE(metallicity, "0.04]")
        ATTRIBUTE_DEFAULT_VALUE(metallicity, "0.02")

        PROPERTY_DOUBLE(SFE, "star formation efficiency ")
        ATTRIBUTE_MIN_VALUE(SFE, "[.01")
        ATTRIBUTE_MAX_VALUE(SFE, ".15]")
        ATTRIBUTE_DEFAULT_VALUE(SFE, ".025")

        PROPERTY_DOUBLE(cloudNumDensity, "initial cloud number density")
        ATTRIBUTE_QUANTITY(cloudNumDensity, "numbervolumedensity")
        ATTRIBUTE_MIN_VALUE(cloudNumDensity, "[10 /cm3")
        ATTRIBUTE_MAX_VALUE(cloudNumDensity, "2560 /cm3]")
        ATTRIBUTE_DEFAULT_VALUE(cloudNumDensity, "320 /cm3")

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
