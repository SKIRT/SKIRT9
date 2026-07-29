/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TIMEINSTRUMENT_HPP
#define TIMEINSTRUMENT_HPP

#include "ApertureInstrument.hpp"
#include "TimeGrid.hpp"

////////////////////////////////////////////////////////////////////

/** TimeInstrument is an abstract subclass of ApertureInstrument (and thus of DistantInstrument)
    that offers a TimeGrid property, used by subclasses that record the timelag response to a pulse
    in the source luminosity. The origin of the arrival time is defined as the time at which a
    photon emitted at the spatial origin arrives at the observer directly. */
class TimeInstrument : public ApertureInstrument
{
    ITEM_ABSTRACT(TimeInstrument, ApertureInstrument,
                  "an instrument that records the timelag response to a source luminosity pulse")
        ATTRIBUTE_TYPE_ALLOWED_IF(TimeInstrument, "ExtinctionOnly")
        ATTRIBUTE_TYPE_DISPLAYED_IF(TimeInstrument, "Level2")

        PROPERTY_ITEM(timeGrid, TimeGrid, "the time grid for this instrument")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function configures the FluxRecorder instance associated with this instrument. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
