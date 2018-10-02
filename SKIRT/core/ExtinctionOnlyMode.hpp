/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef EXTINCTIONONLYMODE_HPP
#define EXTINCTIONONLYMODE_HPP

#include "WithMediumMode.hpp"

////////////////////////////////////////////////////////////////////

/** ExtinctionOnlyMode is a SimulationMode subclass indicating a simulation mode that calculates
    the extinction of the primary radiation through the configured media, including the effects of
    absorption and scattering. In this mode, the simulation does not track the radation field nor
    the state of any media properties. This simulation mode is meaningful only for wavelengths at
    which secondary sources (radiation from the media) can be neglected, i.e. in the ultraviolet,
    optical and near-infrared. */
class ExtinctionOnlyMode : public WithMediumMode
{
    ITEM_CONCRETE(ExtinctionOnlyMode, WithMediumMode, "an extinction-only simulation (no secondary emission)")
    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif
