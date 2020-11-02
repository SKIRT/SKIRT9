/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef EXTINCTIONONLYOPTIONS_HPP
#define EXTINCTIONONLYOPTIONS_HPP

#include "DisjointWavelengthGrid.hpp"
#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** The ExtinctionOnlyOptions class simply offers a number of configuration options that are
    relevant for one of the "extinction only" simulation modes. In these modes, there is no need to
    store the radiation field during the photon packet life cycle. However, there is
    user-configurable option to store the radiation field anyway so that it can be probed for
    output. */
class ExtinctionOnlyOptions : public SimulationItem
{
    ITEM_CONCRETE(ExtinctionOnlyOptions, SimulationItem, "a set of options related to extinction-only simulation modes")

        PROPERTY_BOOL(storeRadiationField, "store the radiation field so that it can be probed for output")
        ATTRIBUTE_DEFAULT_VALUE(storeRadiationField, "false")
        ATTRIBUTE_RELEVANT_IF(storeRadiationField, "ForceScattering")
        ATTRIBUTE_DISPLAYED_IF(storeRadiationField, "Level3")
        ATTRIBUTE_INSERT(storeRadiationField, "storeRadiationField:RadiationField")

        PROPERTY_ITEM(radiationFieldWLG, DisjointWavelengthGrid, "the wavelength grid for storing the radiation field")
        ATTRIBUTE_RELEVANT_IF(radiationFieldWLG, "ForceScattering&storeRadiationField&Panchromatic")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif
