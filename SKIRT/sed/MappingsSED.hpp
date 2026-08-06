/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MAPPINGSSED_HPP
#define MAPPINGSSED_HPP

#include "FamilySED.hpp"

////////////////////////////////////////////////////////////////////

/** MappingsSED is a class that represents spectral energy distributions of star-forming regions
    according to the Mappings III model, parameterized on metallicity, compactness, ISM pressure
    and PDR covering factor. See the MappingsSEDFamily class for more information. */
class MappingsSED : public FamilySED
{
    ITEM_CONCRETE(MappingsSED, FamilySED, "a star-forming region SED from the MAPPINGS III model")
        ATTRIBUTE_TYPE_DISPLAYED_IF(MappingsSED, "Level2")

        PROPERTY_DOUBLE(metallicity, "the metallicity")
        ATTRIBUTE_MIN_VALUE(metallicity, "[0.0006")
        ATTRIBUTE_MAX_VALUE(metallicity, "0.025]")
        ATTRIBUTE_DEFAULT_VALUE(metallicity, "0.0122")

        PROPERTY_DOUBLE(compactness, "the logarithm of the compactness parameter")
        ATTRIBUTE_MIN_VALUE(compactness, "[4.0")
        ATTRIBUTE_MAX_VALUE(compactness, "6.5]")
        ATTRIBUTE_DEFAULT_VALUE(compactness, "6.0")

        PROPERTY_DOUBLE(pressure, "the ISM pressure")
        ATTRIBUTE_QUANTITY(pressure, "pressure")
        ATTRIBUTE_MIN_VALUE(pressure, "[1e10 K/m3")
        ATTRIBUTE_MAX_VALUE(pressure, "1e14 K/m3]")
        ATTRIBUTE_DEFAULT_VALUE(pressure, "1e11 K/m3")

        PROPERTY_DOUBLE(coveringFactor, "the PDR covering factor")
        ATTRIBUTE_MIN_VALUE(coveringFactor, "[0")
        ATTRIBUTE_MAX_VALUE(coveringFactor, "1]")
        ATTRIBUTE_DEFAULT_VALUE(coveringFactor, "0.2")

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
