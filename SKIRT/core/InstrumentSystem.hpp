/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef INSTRUMENTSYSTEM_HPP
#define INSTRUMENTSYSTEM_HPP

#include "Instrument.hpp"
#include "WavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** An InstrumentSystem instance keeps a list of zero or more instruments and an optional default
    wavelength grid that will be used by an instrument unless it specifies its own wavelength grid.
    The instruments can be of various nature and do not need to be located at the same observing
    position. */
class InstrumentSystem : public SimulationItem
{
    ITEM_CONCRETE(InstrumentSystem, SimulationItem, "an instrument system")

        PROPERTY_ITEM(defaultWavelengthGrid, WavelengthGrid, "the default instrument wavelength grid")
        ATTRIBUTE_DEFAULT_VALUE(defaultWavelengthGrid, "LogWavelengthGrid")
        ATTRIBUTE_RELEVANT_IF(defaultWavelengthGrid, "Panchromatic")
        ATTRIBUTE_REQUIRED_IF(defaultWavelengthGrid, "!Level2")
        ATTRIBUTE_INSERT(defaultWavelengthGrid, "defaultWavelengthGrid:DefaultInstrumentWavelengthGrid")

        PROPERTY_ITEM_LIST(instruments, Instrument, "the instruments")
        ATTRIBUTE_DEFAULT_VALUE(instruments, "SEDInstrument")
        ATTRIBUTE_REQUIRED_IF(instruments, "false")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function calls the determineSameObserverAsPreceding() function for all instruments in
        the instrument system except for the first one (because it doesn't have a preceding
        instrument). */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function flushes any information buffered during photon packet detection for the
        complete instrument system. It calls the flush() function for each of the instruments. */
    void flush();

    /** This function writes the recorded data for the complete instrument system to a set of
        files. It calls the write() function for each of the instruments. */
    void write();
};

////////////////////////////////////////////////////////////////////

#endif
