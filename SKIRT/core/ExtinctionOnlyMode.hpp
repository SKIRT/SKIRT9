/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef EXTINCTIONONLYMODE_HPP
#define EXTINCTIONONLYMODE_HPP

#include "WithMediumMode.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

/** ExtinctionOnlyMode is a SimulationMode subclass indicating a simulation mode that calculates
    the extinction of the primary radiation through the configured media, including the effects of
    absorption and scattering. This simulation mode is meaningful only for wavelengths at which
    secondary sources (radiation from the media) can be neglected, i.e. in the ultraviolet, optical
    and near-infrared.

    In this simulation mode, there is no secondary emission, and the media state is constant, i.e.
    it is initialized from the properties defined in the input model (density distributions,
    material properties) and is never updated. As a result, there is no need to store the radiation
    field during the photon packet life cycle. However, there is user-configurable option to store
    the radiation field anyway so that it can be probed for output. */
class ExtinctionOnlyMode : public WithMediumMode
{
    ITEM_CONCRETE(ExtinctionOnlyMode, WithMediumMode, "an extinction-only simulation (no secondary emission)")
        ATTRIBUTE_TYPE_INSERT(ExtinctionOnlyMode, "ExtinctionOnly")

    PROPERTY_BOOL(storeRadiationField, "store the radiation field so that it can be probed for output")
        ATTRIBUTE_DEFAULT_VALUE(storeRadiationField, "false")
        ATTRIBUTE_DISPLAYED_IF(storeRadiationField, "Level3")
        ATTRIBUTE_INSERT(storeRadiationField, "storeRadiationField:RadiationField")

    PROPERTY_ITEM(radiationFieldWLG, WavelengthGrid, "the wavelength grid for storing the radiation field")
        ATTRIBUTE_RELEVANT_IF(radiationFieldWLG, "storeRadiationField&Panchromatic")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif
