/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Instrument.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "FluxRecorder.hpp"

////////////////////////////////////////////////////////////////////

void Instrument::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // select "local" or default wavelength grid
    _instrumentWavelengthGrid = wavelengthGrid() ? wavelengthGrid() : find<WavelengthGrid>();

    // TO DO: discover details about the simulation
    auto config = find<Configuration>();
    bool hasMedium = config->hasMedium();
    bool hasMediumEmission = false;

    // partially configure the flux recorder
    _recorder = new FluxRecorder(this);
    _recorder->setSimulationInfo(instrumentName(), instrumentWavelengthGrid(), hasMedium, hasMediumEmission);
    _recorder->setUserFlags(_recordComponents, _numScatteringLevels, _recordPolarization, _recordStatistics);
}

////////////////////////////////////////////////////////////////////

void Instrument::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    // finalize configuration of the flux recorder
    _recorder->finalizeConfiguration();
}

////////////////////////////////////////////////////////////////////

Instrument::~Instrument()
{
    delete _recorder;
}

////////////////////////////////////////////////////////////////////

std::string Instrument::itemName() const
{
    return instrumentName();
}

////////////////////////////////////////////////////////////////////

void Instrument::flush()
{
    _recorder->flush();
}

////////////////////////////////////////////////////////////////////

void Instrument::write()
{
    _recorder->calibrateAndWrite();
}

////////////////////////////////////////////////////////////////////
