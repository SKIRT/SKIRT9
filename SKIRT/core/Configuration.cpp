/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Configuration.hpp"
#include "DustEmissionMode.hpp"
#include "ExtinctionOnlyMode.hpp"
#include "FatalError.hpp"
#include "MaterialMix.hpp"
#include "MonteCarloSimulation.hpp"
#include "NR.hpp"
#include "NoMediumMode.hpp"
#include "OligoWavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

Configuration::Configuration(SimulationItem* parent)
{
    parent->addChild(this);
}

////////////////////////////////////////////////////////////////////

void Configuration::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // Implementation note: this function is NOT allowed to perform setup on any simulation item in the hierarchy;
    //                      in other words, always use find<XXX>(false) and check for a nullptr result.

    // retrieve objects that we'll need anyway
    auto sim = find<MonteCarloSimulation>(false);
    if (!sim) throw FATALERROR("Cannot locate a MonteCarloSimulation object in the simulation hierarchy");
    auto mode = find<SimulationMode>(false);
    if (!mode) throw FATALERROR("Cannot locate a SimulationMode object in the simulation hierarchy");
    auto ss = find<SourceSystem>(false);
    if (!ss) throw FATALERROR("Cannot locate a SourceSystem object in the simulation hierarchy");

    // retrieve wavelength-related options
    _oligochromatic = sim->wavelengthRegime() == MonteCarloSimulation::WavelengthRegime::Oligochromatic;
    if (_oligochromatic)
    {
        NR::assign(_oligoWavelengths, ss->wavelengths());
        _defaultWavelengthGrid = new OligoWavelengthGrid(this, ss->wavelengths());
        _oligoBinWidth = _defaultWavelengthGrid->effectiveWidth(0);
        _sourceWavelengthRange.set(_defaultWavelengthGrid->leftBorder(0),
                                   _defaultWavelengthGrid->rightBorder(_defaultWavelengthGrid->numBins()-1));
    }
    else
    {
        auto is = find<InstrumentSystem>(false);
        if (is) _defaultWavelengthGrid = is->defaultWavelengthGrid();
        _sourceWavelengthRange.set(ss->minWavelength(), ss->maxWavelength());
    }

    // retrieve photon life-cycle and medium-related options
    _numPrimaryPackets = sim->numPackets();
    bool mustHaveMedium = false;
    auto mmode = dynamic_cast<WithMediumMode*>(mode);
    if (mmode)
    {
        mustHaveMedium = true;
        _minWeightReduction = mmode->minWeightReduction();
        _minScattEvents = mmode->minScattEvents();
        _pathLengthBias = mmode->pathLengthBias();
        _numDensitySamples = mmode->numDensitySamples();
    }

    // retrieve dust emission-related options
    auto emode = dynamic_cast<DustEmissionMode*>(mode);
    if (emode)
    {
        _numPrimaryPackets = sim->numPackets() * emode->primaryPacketsMultiplier();
        _numSecondaryPackets = sim->numPackets() * emode->secondaryPacketsMultiplier();
        throw FATALERROR("Dust emission is not yet supported");
    }

    // determine the number of media in the simulation hierarchy
    int numMedia = 0;
    auto ms = find<MediumSystem>(false);
    if (ms) numMedia = ms->media().size();  // may be zero
    _hasMedium = (numMedia!=0);

    // verify this with the requirements set by the simulation mode
    if (!mustHaveMedium && _hasMedium)
        throw FATALERROR("This simulation mode does not allow media to be configured");
    if (mustHaveMedium && !_hasMedium)
        throw FATALERROR("This simulation mode requires at least one medium to be configured");

    // check for polarization
    if (_hasMedium)
    {
        int numPolarization = 0;
        for (auto medium : ms->media()) if (medium->mix()->hasPolarization()) numPolarization++;
        if (numPolarization!=0 && numPolarization!=numMedia)
            throw FATALERROR("All media must consistenly support polarization, or not support polarization");
        _hasPolarization = numPolarization!=0;
    }

    // TO DO: retrieve options from and verify configuration for other modes

    // in case emulation mode has been set before our setup() was called, perform the emulation overrides again
    if (emulationMode()) setEmulationMode();
}

////////////////////////////////////////////////////////////////////

void Configuration::setEmulationMode()
{
    _emulationMode = true;
    _numPrimaryPackets = 0.;
    // TO DO: set other number of packet variables to zero
    // TO DO: force the number of state iterations to one, if applicable
}

////////////////////////////////////////////////////////////////////

WavelengthGrid* Configuration::wavelengthGrid(WavelengthGrid* localWavelengthGrid) const
{
    auto result = localWavelengthGrid && !_oligochromatic ? localWavelengthGrid : _defaultWavelengthGrid;
    if (!result) throw FATALERROR("Cannot find a wavelength grid for instrument or probe");
    return result;
}

////////////////////////////////////////////////////////////////////
