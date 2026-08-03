/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SELECTDUSTMIXFAMILY_HPP
#define SELECTDUSTMIXFAMILY_HPP

#include "DustMix.hpp"
#include "MaterialMixFamily.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the SelectDustMixFamily class represents a family of dust mixes that is
    specified as part of the configuration. Specifically, a property of this class holds a
    user-configurable list of dust mixes representing the family. The family requires a single
    parameter value to select a family member, corresponding to the zero-based index in the
    configured list of dust mixes. The floating point parameter value is rounded to the nearest
    integer and subsequently clipped to be in range. */
class SelectDustMixFamily : public MaterialMixFamily
{
    ITEM_CONCRETE(SelectDustMixFamily, MaterialMixFamily, "a family of dust mixes specified in the configuration")

        PROPERTY_ITEM_LIST(dustMixes, DustMix, "the family of dust mixes")
        ATTRIBUTE_DEFAULT_VALUE(dustMixes, "ConfigurableDustMix")

    ITEM_END()

    //====================== Other functions =====================

public:
    /** This function returns SnapshotParameter information for a single parameter corresponding to
        the zero-based index in the configured list of dust mixes. */
    vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns (a pointer to) the dust mix corresponding to the zero-based index in
        the configured list of dust mixes specified as the single parameter in the list. The
        floating point parameter value is rounded to the nearest integer and subsequently clipped
        to be in range. If the number of parameters in the specified list is not equal to one, the
        behavior is undefined. The material mix family retains ownership of the returned dust mix,
        and guarantees that it will not be destroyed until the family itself is destroyed. */
    const MaterialMix* mix(const Array& parameters) const override;
};

////////////////////////////////////////////////////////////////////

#endif
