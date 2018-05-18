/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef INSTRUMENT_HPP
#define INSTRUMENT_HPP

#include "SimulationItem.hpp"
#include "Direction.hpp"
#include "FluxRecorder.hpp"
#include "Position.hpp"
#include "WavelengthGrid.hpp"
class PhotonPacket;

////////////////////////////////////////////////////////////////////

/** Instrument is an abstract class representing instruments to collect the photon packets during a
    Monte Carlo simulation. Subclasses implement instruments that can vary in type of projection
    (e.g. parallel or perspective) and in what is recorded (e.g. SED or full integral field data
    cube). This top-level abstract class offers a generic interface for receiving photon packets
    from the simulation. It also includes facilities for configuring user properties that are
    common to all instruments, such as which flux contributions need to be recorded. A wavelength
    grid is established either by specifying a grid for this instrument specifically, or by
    defaulting to the common grid specified for the instrument system. */
class Instrument : public SimulationItem
{
    ITEM_ABSTRACT(Instrument, SimulationItem, "an instrument")

    PROPERTY_STRING(instrumentName, "the name for this instrument")

    PROPERTY_ITEM(wavelengthGrid, WavelengthGrid, "the wavelength grid for this instrument")
        ATTRIBUTE_OPTIONAL(wavelengthGrid)

    PROPERTY_BOOL(recordComponents, "record flux components separately")
        ATTRIBUTE_DEFAULT_VALUE(recordComponents, "false")

    PROPERTY_INT(numScatteringLevels, "the number of individually recorded scattering levels")
        ATTRIBUTE_RELEVANT_IF(numScatteringLevels, "recordComponents")
        ATTRIBUTE_MIN_VALUE(numScatteringLevels, "0")
        ATTRIBUTE_MAX_VALUE(numScatteringLevels, "99")
        ATTRIBUTE_DEFAULT_VALUE(numScatteringLevels, "0")

    PROPERTY_BOOL(recordPolarization, "record polarization (Stokes vector elements)")
        ATTRIBUTE_DEFAULT_VALUE(recordPolarization, "false")

    PROPERTY_BOOL(recordStatistics, "record information for calculating statistical properties")
        ATTRIBUTE_DEFAULT_VALUE(recordStatistics, "false")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function performs setup for the instrument. It establishes the wavelength grid for the
        instrument: if a grid is specified for the instrument, that grid is used. If not, the
        default wavelength grid specified for the instrument system is used instead. If neither of
        these grids are specified, the function throws a fatal error.

        The function also partially configures the FluxRecorder instance for this instrument,
        passing it the values of the user properties offered by this class and some extra
        information on the simulation. The setupSelfBefore() function of each subclass is expected
        to augment the configuration by calling the includeFluxDensity() and/or
        includeSurfaceBrightness() functions. */
    void setupSelfBefore() override;

    /** This function finalizes the configuration of FluxRecorder instance for this instrument. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** Returns the wavelength grid for the instrument as determined during setup, i.e. either the
        grid specified for this instrument or the default grid specified for the instrument system.
        After setup has completed, the function never returns a nulltpr because setupSelfBefore()
        throws a fatal error if neither of these grids are specified. */
    const WavelengthGrid* instrumentWavelengthGrid() const { return _instrumentWavelengthGrid; }

    /** This function flushes any information buffered by the detect() function. It simply calls
        the corresponding function of the FluxRecorder instance associated with this instrument. */
    void flush();

    /** This function calibrates the instrument and outputs the recorded contents to a set of
        files. It simply calls the corresponding function of the FluxRecorder instance associated
        with this instrument. */
    void write();

protected:
    /** Returns the FluxRecorder instance associated with this instrument. This function is
        intended for use in subclasses only. */
    FluxRecorder* instrumentFluxRecorder() { return &_recorder; }

    //=========== Functions to be implemented in subclass ===========

public:
    /** Returns the direction towards the observer, given the photon packet's launching position.
        The implementation must be provided in a subclass. */
    virtual Direction bfkobs(const Position& bfr) const = 0;

    /** Returns the direction along the positive x-axis of the instrument frame, expressed in model
        coordinates. The implementation must be provided in a subclass. */
    virtual Direction bfkx() const = 0;

    /** Returns the direction along the positive y-axis of the instrument frame, expressed in model
        coordinates. The implementation must be provided in a subclass. */
    virtual Direction bfky() const = 0;

    /** This function simulates the detection of a photon packet by the instrument. Its
        implementation must be provided in a subclass. */
    virtual void detect(PhotonPacket* pp) = 0;

    //======================== Data Members =======================

private:
    const WavelengthGrid* _instrumentWavelengthGrid{nullptr};
    FluxRecorder _recorder;
};

////////////////////////////////////////////////////////////////////

#endif
