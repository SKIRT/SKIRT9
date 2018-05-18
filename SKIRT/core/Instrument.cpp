/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Instrument.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

void Instrument::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // select "local" or default wavelength grid
    _instrumentWavelengthGrid = wavelengthGrid() ? wavelengthGrid() : find<WavelengthGrid>();

    // TO DO: discover details about the simulation
    bool hasMedium = false;
    bool hasMediumEmission = false;

    // partially configure the flux recorder
    _recorder.setSimulationInfo(instrumentName(), instrumentWavelengthGrid(), hasMedium, hasMediumEmission);
    _recorder.setUserFlags(_recordComponents, _numScatteringLevels, _recordPolarization, _recordStatistics);
}

////////////////////////////////////////////////////////////////////

void Instrument::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    // finalize configuration of the flux recorder
    _recorder.finalizeConfiguration();
}

////////////////////////////////////////////////////////////////////

void Instrument::flush()
{
    _recorder.flush();
}

////////////////////////////////////////////////////////////////////

void Instrument::write()
{
    _recorder.calibrateAndWrite();
}

////////////////////////////////////////////////////////////////////
