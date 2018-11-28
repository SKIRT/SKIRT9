/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Configuration.hpp"
#include "AllCellsLibrary.hpp"
#include "DustEmissionOptions.hpp"
#include "DustSelfAbsorptionOptions.hpp"
#include "ExtinctionOnlyOptions.hpp"
#include "FatalError.hpp"
#include "MaterialMix.hpp"
#include "MonteCarloSimulation.hpp"
#include "NR.hpp"
#include "OligoWavelengthDistribution.hpp"
#include "OligoWavelengthGrid.hpp"
#include "PhotonPacketOptions.hpp"

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
    auto ss = find<SourceSystem>(false);
    if (!ss) throw FATALERROR("Cannot locate a SourceSystem object in the simulation hierarchy");
    auto is = find<InstrumentSystem>(false);

    // retrieve wavelength-related options
    _oligochromatic = sim->simulationMode() == MonteCarloSimulation::SimulationMode::OligoNoMedium ||
                      sim->simulationMode() == MonteCarloSimulation::SimulationMode::OligoExtinctionOnly;
    if (_oligochromatic)
    {
        auto oligoWavelengthGrid = new OligoWavelengthGrid(this, ss->wavelengths());
        _sourceWavelengthRange = oligoWavelengthGrid->wavelengthRange();
        _defaultWavelengthGrid = oligoWavelengthGrid;
        _oligoWavelengthBiasDistribution = new OligoWavelengthDistribution(oligoWavelengthGrid);
    }
    else
    {
        _sourceWavelengthRange.set(ss->minWavelength(), ss->maxWavelength());
        if (is) _defaultWavelengthGrid = is->defaultWavelengthGrid();
    }

    // determine the number of media in the simulation hierarchy
    int numMedia = 0;
    auto ms = find<MediumSystem>(false);
    if (ms) numMedia = ms->media().size();  // may be zero
    _hasMedium = (numMedia!=0);

    // verify this with the requirements set by the simulation mode
    bool mustHaveMedium = sim->simulationMode() != MonteCarloSimulation::SimulationMode::OligoNoMedium &&
                          sim->simulationMode() != MonteCarloSimulation::SimulationMode::NoMedium;;
    if (!mustHaveMedium && _hasMedium)
        throw FATALERROR("This simulation mode does not allow media to be configured");
    if (mustHaveMedium && !_hasMedium)
        throw FATALERROR("This simulation mode requires at least one medium to be configured");

    // retrieve photon life-cycle and basic medium-related options
    _numPrimaryPackets = sim->numPackets();
    if (_hasMedium)
    {
        _numDensitySamples = ms->numDensitySamples();
        _minWeightReduction = ms->photonPacketOptions()->minWeightReduction();
        _minScattEvents = ms->photonPacketOptions()->minScattEvents();
        _pathLengthBias = ms->photonPacketOptions()->pathLengthBias();
    }

    // retrieve extinction-only options
    if (sim->simulationMode() == MonteCarloSimulation::SimulationMode::OligoExtinctionOnly ||
        sim->simulationMode() == MonteCarloSimulation::SimulationMode::ExtinctionOnly)
    {
        _hasRadiationField = ms->extinctionOnlyOptions()->storeRadiationField();
        if (_hasRadiationField)
            _radiationFieldWLG = _oligochromatic ? dynamic_cast<OligoWavelengthGrid*>(_defaultWavelengthGrid)
                                                 : ms->extinctionOnlyOptions()->radiationFieldWLG();
    }

    // retrieve dust emission options
    if (sim->simulationMode() == MonteCarloSimulation::SimulationMode::DustEmission ||
        sim->simulationMode() == MonteCarloSimulation::SimulationMode::DustEmissionWithSelfAbsorption)
    {
        _hasRadiationField = true;
        _hasDustEmission = true;
        _dustEmissivity = ms->dustEmissionOptions()->dustEmissivity();
        _cellLibrary = ms->dustEmissionOptions()->cellLibrary();
        if (!_cellLibrary) _cellLibrary = new AllCellsLibrary(this);
        _radiationFieldWLG = ms->dustEmissionOptions()->radiationFieldWLG();
        _dustEmissionWLG = ms->dustEmissionOptions()->dustEmissionWLG();
        _numSecondaryPackets = sim->numPackets() * ms->dustEmissionOptions()->secondaryPacketsMultiplier();
        _secondarySpatialBias = ms->dustEmissionOptions()->spatialBias();
        _secondaryWavelengthBias = ms->dustEmissionOptions()->wavelengthBias();
        _secondaryWavelengthBiasDistribution = ms->dustEmissionOptions()->wavelengthBiasDistribution();
    }

    // retrieve dust self-absorption options
    if (sim->simulationMode() == MonteCarloSimulation::SimulationMode::DustEmissionWithSelfAbsorption)
    {
        _hasDustSelfAbsorption = true;
        _minIterations = ms->dustSelfAbsorptionOptions()->minIterations();
        _maxIterations = ms->dustSelfAbsorptionOptions()->maxIterations();
        _maxFractionOfPrimary = ms->dustSelfAbsorptionOptions()->maxFractionOfPrimary();
        _maxFractionOfPrevious = ms->dustSelfAbsorptionOptions()->maxFractionOfPrevious();
        _numIterationPackets = sim->numPackets() * ms->dustSelfAbsorptionOptions()->iterationPacketsMultiplier();
    }

    // retrieve symmetry dimensions
    if (_hasMedium)
    {
        _modelDimension = max(ss->dimension(), ms->dimension());
        _gridDimension = ms->gridDimension();
        if (_modelDimension > _gridDimension)
            throw FATALERROR("The grid symmetry (" + std::to_string(_gridDimension) + "D)"
                             "does not support the model symmetry (" + std::to_string(_modelDimension) + "D)");
    }
    else
    {
        _modelDimension = ss->dimension();
    }

    // check for polarization
    if (_hasMedium)
    {
        int numPolarization = 0;
        for (auto medium : ms->media()) if (medium->mix()->hasPolarization()) numPolarization++;
        if (numPolarization!=0 && numPolarization!=numMedia)
            throw FATALERROR("All media must consistenly support polarization, or not support polarization");
        _hasPolarization = numPolarization!=0;
    }

    // check for velocities in media
    if (_hasMedium) for (auto medium : ms->media()) if (medium->hasVelocity()) _hasMovingMedia = true;

    // check for variable material mixes
    // TO DO: implement this check once variable material mixes have been actually implemented
    // TO DO: warn when using cell libraries in combination with variable material mixes

    // in case emulation mode has been set before our setup() was called, perform the emulation overrides again
    if (emulationMode()) setEmulationMode();
}

////////////////////////////////////////////////////////////////////

void Configuration::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    // Implementation note: this function should perform pure logging and is NOT allowed to perform setup
    //                      on any simulation item in the hierarchy except for the logger.
    auto log = find<Log>();

    // --- log model symmetries ---

    // if there are no media, simply log the source model symmetry
    if (!_hasMedium)
    {
        log->info("Model symmetry: " + std::to_string(_modelDimension) + "D");
    }

    // if there are media, compare the model symmetry to the grid symmetry
    // (the case where the grid has insufficient dimension causes a fatal error in setupSelfBefore)
    else
    {
        if (_modelDimension == _gridDimension)
        {
            log->info("Model and grid symmetry: " + std::to_string(_modelDimension) + "D");
        }
        else
        {
            log->info("Model symmetry: " + std::to_string(_modelDimension) + "D; "
                      "Spatial grid symmetry: " + std::to_string(_gridDimension) + "D");
            log->warning("Selecting a grid with the model symmetry might be more efficient");
        }
    }
}

////////////////////////////////////////////////////////////////////

void Configuration::setEmulationMode()
{
    _emulationMode = true;
    _numPrimaryPackets = 0.;
    _numIterationPackets = 0.;
    _numSecondaryPackets = 0.;
    _minIterations = 1;
    _maxIterations = 1;
}

////////////////////////////////////////////////////////////////////

Range Configuration::simulationWavelengthRange() const
{
    // include primary and secondary source ranges
    Range range = _sourceWavelengthRange;
    if (_dustEmissionWLG)
    {
        _dustEmissionWLG->setup();      // ensure setup because this function may be called early during setup
        range.extend(_dustEmissionWLG->wavelengthRange());
    }

    // extend the range with a wide margin for kinematics
    double z = 1./3.;
    return Range( range.min()/(1.+z), range.max()*(1.+z) );
}

////////////////////////////////////////////////////////////////////

WavelengthGrid* Configuration::wavelengthGrid(WavelengthGrid* localWavelengthGrid) const
{
    auto result = localWavelengthGrid && !_oligochromatic ? localWavelengthGrid : _defaultWavelengthGrid;
    if (!result) throw FATALERROR("Cannot find a wavelength grid for instrument or probe");
    return result;
}

////////////////////////////////////////////////////////////////////
