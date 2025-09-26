/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Instrument.hpp"
#include "Configuration.hpp"
#include "FluxRecorder.hpp"

////////////////////////////////////////////////////////////////////

void Instrument::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // select "local" or default wavelength grid
    auto config = find<Configuration>();
    _instrumentWavelengthGrid = config->wavelengthGrid(wavelengthGrid());

    // discover details about the simulation
    bool hasMedium = config->hasMedium();
    bool hasMediumEmission = config->hasSecondaryEmission() || config->scatteringEmulatesSecondaryEmission();

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
