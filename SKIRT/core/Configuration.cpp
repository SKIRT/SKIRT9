/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Configuration.hpp"
#include "AllCellsLibrary.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"
#include "MaterialMix.hpp"
#include "MaterialWavelengthRangeInterface.hpp"
#include "MonteCarloSimulation.hpp"
#include "NR.hpp"
#include "OligoWavelengthDistribution.hpp"
#include "OligoWavelengthGrid.hpp"
#include "StringUtils.hpp"
#include "VoronoiMeshSpatialGrid.hpp"
#include <set>

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

    // locate objects that we'll need anyway
    auto sim = find<MonteCarloSimulation>(false);
    if (!sim) throw FATALERROR("Cannot locate a MonteCarloSimulation object in the simulation hierarchy");
    auto ss = find<SourceSystem>(false);
    if (!ss) throw FATALERROR("Cannot locate a SourceSystem object in the simulation hierarchy");

    // ---- set configuration relevant for all simulations, even if there are no media ----

    // retrieve simulation mode
    auto simulationMode = sim->simulationMode();
    using SimulationMode = MonteCarloSimulation::SimulationMode;

    // retrieve base number of packets
    _numPrimaryPackets = sim->numPackets();

    // retrieve cosmology parameters
    _redshift = sim->cosmology()->modelRedshift();
    _angularDiameterDistance = sim->cosmology()->angularDiameterDistance();
    _luminosityDistance = sim->cosmology()->luminosityDistance();

    // retrieve wavelength-related options
    _oligochromatic =
        simulationMode == SimulationMode::OligoNoMedium || simulationMode == SimulationMode::OligoExtinctionOnly;
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
        auto is = find<InstrumentSystem>(false);
        if (is) _defaultWavelengthGrid = is->defaultWavelengthGrid();
    }

    // determine model dimension based on sources only (we'll redo this later if there are media)
    _modelDimension = ss->dimension();

    // determine the number of media in the simulation hierarchy
    int numMedia = 0;
    auto ms = find<MediumSystem>(false);
    if (ms) numMedia = ms->media().size();  // may be zero
    _hasMedium = (numMedia != 0);

    // verify this with the requirements set by the simulation mode
    bool mustHaveMedium = simulationMode != SimulationMode::OligoNoMedium && simulationMode != SimulationMode::NoMedium;
    if (!mustHaveMedium && _hasMedium) throw FATALERROR("This simulation mode does not allow media to be configured");
    if (mustHaveMedium && !_hasMedium)
        throw FATALERROR("This simulation mode requires at least one medium to be configured");

    // if there are no media, we're done
    if (!_hasMedium) return;

    // ---- set configuration relevant for simulations that have media ----

    // retrieve basic photon life-cycle options
    _forceScattering = ms->photonPacketOptions()->forceScattering();
    _minWeightReduction = ms->photonPacketOptions()->minWeightReduction();
    _minScattEvents = ms->photonPacketOptions()->minScattEvents();
    _pathLengthBias = ms->photonPacketOptions()->pathLengthBias();

    // retrieve Lyman-alpha options
    if (simulationMode == SimulationMode::LyaExtinctionOnly)
    {
        _hasLymanAlpha = true;
        switch (ms->lyaOptions()->lyaAccelerationScheme())
        {
            case LyaOptions::LyaAccelerationScheme::None:
                _lyaAccelerationScheme = Configuration::LyaAccelerationScheme::None;
                break;
            case LyaOptions::LyaAccelerationScheme::Constant:
                _lyaAccelerationScheme = Configuration::LyaAccelerationScheme::Constant;
                _lyaAccelerationStrength = ms->lyaOptions()->lyaAccelerationStrength();
                break;
            case LyaOptions::LyaAccelerationScheme::Variable:
                _lyaAccelerationScheme = Configuration::LyaAccelerationScheme::Variable;
                _lyaAccelerationStrength = ms->lyaOptions()->lyaAccelerationStrength();
                break;
        }
        if (ms->lyaOptions()->includeHubbleFlow()) _hubbleExpansionRate = sim->cosmology()->relativeExpansionRate();
    }

    // retrieve the presence of phases and iterations
    if (simulationMode == SimulationMode::ExtinctionOnly)
    {
        _hasPrimaryIterations = sim->iterateMediumState();
    }
    else if (simulationMode == SimulationMode::DustEmission)
    {
        _hasSecondaryEmission = true;
        _hasDustEmission = true;
        _hasPrimaryIterations = sim->iterateMediumState();
        _hasSecondaryIterations = sim->iterateSecondaryEmission();
    }
    else if (simulationMode == SimulationMode::GasEmission)
    {
        _hasSecondaryEmission = true;
        _hasGasEmission = true;
        _hasPrimaryIterations = sim->iterateMediumState();
        _hasSecondaryIterations = sim->iterateSecondaryEmission();
    }
    else if (simulationMode == SimulationMode::DustAndGasEmission)
    {
        _hasSecondaryEmission = true;
        _hasDustEmission = true;
        _hasGasEmission = true;
        _hasPrimaryIterations = sim->iterateMediumState();
        _hasSecondaryIterations = sim->iterateSecondaryEmission();
    }

    // retrieve secondary emission options
    if (_hasSecondaryEmission)
    {
        _numSecondaryPackets = _numPrimaryPackets * ms->secondaryEmissionOptions()->secondaryPacketsMultiplier();
        _storeEmissionRadiationField = ms->secondaryEmissionOptions()->storeEmissionRadiationField();
        _secondarySpatialBias = ms->secondaryEmissionOptions()->spatialBias();
        _secondarySourceBias = ms->secondaryEmissionOptions()->sourceBias();
    }

    // retrieve primary iteration options; suppress primary iteration if no packets are launched
    if (_hasPrimaryIterations)
    {
        _numPrimaryIterationPackets = _numPrimaryPackets * ms->iterationOptions()->iterationPacketsMultiplier();
        if (_numPrimaryIterationPackets >= 1.)
        {
            _minPrimaryIterations = ms->iterationOptions()->minPrimaryIterations();
            _maxPrimaryIterations = max(_minPrimaryIterations, ms->iterationOptions()->maxPrimaryIterations());
        }
        else
            _hasPrimaryIterations = false;
    }

    // if iteration over primary emission is requested, verify that there are dynamic state recipes
    if (_hasPrimaryIterations)
    {
        if (ms->dynamicStateOptions()->recipes().empty())
            throw FATALERROR("At least one dynamic state recipe must be configured when iterateMediumState is true");
        _hasDynamicState = true;
    }

    // retrieve secondary iteration options; suppress secondary iteration if no packets are launched
    if (_hasSecondaryIterations)
    {
        _numSecondaryIterationPackets = _numSecondaryPackets * ms->iterationOptions()->iterationPacketsMultiplier();
        if (_numSecondaryIterationPackets > 0.)
        {
            _minSecondaryIterations = ms->iterationOptions()->minSecondaryIterations();
            _maxSecondaryIterations = max(_minSecondaryIterations, ms->iterationOptions()->maxSecondaryIterations());
            if (_hasPrimaryIterations) _includePrimaryEmission = ms->iterationOptions()->includePrimaryEmission();
        }
        else
            _hasSecondaryIterations = false;
    }

    // retrieve radiation field options
    _hasRadiationField =
        _hasPrimaryIterations || _hasSecondaryEmission || ms->radiationFieldOptions()->storeRadiationField();
    if (_hasRadiationField)
    {
        _hasPanRadiationField = !_oligochromatic;
        _radiationFieldWLG = _oligochromatic ? dynamic_cast<OligoWavelengthGrid*>(_defaultWavelengthGrid)
                                             : ms->radiationFieldOptions()->radiationFieldWLG();
    }
    _hasSecondaryRadiationField = _hasSecondaryIterations || _storeEmissionRadiationField;

    // retrieve dust emission options
    if (_hasDustEmission)
    {
        if (ms->dustEmissionOptions()->dustEmissionType() == DustEmissionOptions::EmissionType::Stochastic)
        {
            // verify that all dust mixes are multi-grain when requesting stochastic heating
            for (auto medium : ms->media())
                if (medium->mix()->isDust() && !medium->mix()->hasStochasticDustEmission())
                    throw FATALERROR("When requesting stochastic heating, all dust mixes must be multi-grain");
            _hasStochasticDustEmission = true;
        }
        _includeHeatingByCMB = ms->dustEmissionOptions()->includeHeatingByCMB();
        _cellLibrary = ms->dustEmissionOptions()->cellLibrary();
        if (!_cellLibrary) _cellLibrary = new AllCellsLibrary(this);
        _dustEmissionWLG = ms->dustEmissionOptions()->dustEmissionWLG();
        if (_hasSecondaryIterations)
        {
            _maxFractionOfPrimary = ms->dustEmissionOptions()->maxFractionOfPrimary();
            _maxFractionOfPrevious = ms->dustEmissionOptions()->maxFractionOfPrevious();
        }
        _dustEmissionSourceWeight = ms->dustEmissionOptions()->sourceWeight();
        _dustEmissionWavelengthBias = ms->dustEmissionOptions()->wavelengthBias();
        _dustEmissionWavelengthBiasDistribution = ms->dustEmissionOptions()->wavelengthBiasDistribution();
    }

    // retrieve media sampling options
    _numDensitySamples = ms->samplingOptions()->numDensitySamples();
    _numPropertySamples = ms->samplingOptions()->numPropertySamples();

    // retrieve symmetry dimensions
    _modelDimension = max(ss->dimension(), ms->dimension());
    _gridDimension = ms->gridDimension();
    if (_modelDimension > _gridDimension)
        throw FATALERROR("The grid symmetry (" + std::to_string(_gridDimension)
                         + "D) does not support the model symmetry (" + std::to_string(_modelDimension) + "D)");

    // verify that there is a Lya medium component if required, and none if not required
    int numLyaMedia = 0;
    for (auto medium : ms->media())
        if (medium->mix()->hasResonantScattering()) numLyaMedia++;
    if (_hasLymanAlpha && numLyaMedia < 1)
        throw FATALERROR("Lyman-alpha simulation mode requires a medium component with Lyman-alpha material mix");
    if (!_hasLymanAlpha && numLyaMedia > 0)
        throw FATALERROR("Lyman-alpha material mix is allowed only with Lyman-alpha simulation mode");

    // determine whether media must support the generatePosition() function
    // currently, that function is called only by the VoronoiMeshSpatialGrid class for certain policies
    auto grid = dynamic_cast<VoronoiMeshSpatialGrid*>(ms->grid());
    if (grid)
    {
        auto policy = grid->policy();
        if (policy == VoronoiMeshSpatialGrid::Policy::DustDensity
            || policy == VoronoiMeshSpatialGrid::Policy::ElectronDensity
            || policy == VoronoiMeshSpatialGrid::Policy::GasDensity
            || policy == VoronoiMeshSpatialGrid::Policy::ImportedMesh)
        {
            _mediaNeedGeneratePosition = true;
        }
    }

    // check for semi-dymamic medium state
    if (_hasSecondaryEmission)
        for (auto medium : ms->media())
            if (medium->mix()->hasSemiDynamicMediumState()) _hasSemiDynamicState = true;

    // check for velocities in sources and media
    for (auto source : ss->sources())
        if (source->hasVelocity()) _hasMovingSources = true;
    for (auto medium : ms->media())
        if (medium->hasVelocity()) _hasMovingMedia = true;

    // check for variable material mixes
    for (auto medium : ms->media())
        if (medium->hasVariableMix()) _hasVariableMedia = true;

    // check for dependencies on extra specific state variables
    bool hasExtraSpecificState = false;
    for (auto medium : ms->media())
        if (medium->mix()->hasExtraSpecificState()) hasExtraSpecificState = true;

    // set the combined medium criteria
    _hasConstantPerceivedWavelength = !_hasMovingMedia && !_hubbleExpansionRate;
    bool _hasConstantSectionMedium = _hasConstantPerceivedWavelength && !_hasVariableMedia && !hasExtraSpecificState;
    _hasSingleConstantSectionMedium = numMedia == 1 && _hasConstantSectionMedium;
    _hasMultipleConstantSectionMedia = numMedia > 1 && _hasConstantSectionMedium;

    // check for magnetic fields
    for (auto medium : ms->media())
        if (medium->hasMagneticField()) _hasMagneticField = true;

    // check for polarization
    int numPolarization = 0;
    for (auto medium : ms->media())
        if (medium->mix()->hasPolarizedScattering()) numPolarization++;
    if (numPolarization != 0 && numPolarization != numMedia)
        throw FATALERROR("All media must consistenly support polarization, or not support polarization");
    _hasPolarization = numPolarization != 0;

    // check for polarization by spheroidal particles
    if (_hasPolarization)
        for (auto medium : ms->media())
            if (medium->mix()->hasPolarizedAbsorption() || medium->mix()->hasPolarizedEmission())
                _hasSpheroidalPolarization = true;  // this flag may need to be split over absorption and emission

    // spheroidal particles require a magnetic field
    if (_hasSpheroidalPolarization && !_hasMagneticField)
        throw FATALERROR("Polarization by spheroidal particles requires a magnetic field to determine alignment");

    // prohibit non-identity-mapping cell libraries in combination with spatially varying material mixes
    if ((_hasVariableMedia || hasExtraSpecificState) && _cellLibrary && !dynamic_cast<AllCellsLibrary*>(_cellLibrary))
        throw FATALERROR("Cannot use spatial cell library in combination with spatially varying material mixes");
}

////////////////////////////////////////////////////////////////////

void Configuration::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    // Implementation note: this function should perform pure logging and is NOT allowed to perform setup
    //                      on any simulation item in the hierarchy except for the logger.
    auto log = find<Log>();

    // in case emulation mode has been set before our setup() was called, perform the emulation overrides again
    if (emulationMode()) setEmulationMode();

    // --- log wavelength regime, simulation mode, and media characteristics  ---

    string regime = _oligochromatic ? "Oligo" : "Pan";
    log->info("  " + regime + "chromatic wavelength regime");
    string medium = _hasMedium ? "With" : "No";
    log->info("  " + medium + " transfer medium");
    if (_hasLymanAlpha) log->info("  Including Lyman-alpha line transfer");

    if (_hasPrimaryIterations && _hasSecondaryIterations)
    {
        if (_includePrimaryEmission)
            log->info("  Iterating over primary emission, and then over primary and secondary emission");
        else
            log->info("  Iterating over primary emission, and then over secondary emission");
    }
    else if (_hasPrimaryIterations)
    {
        log->info("  Iterating over primary emission");
    }
    else if (_hasSecondaryIterations)
    {
        log->info("  Iterating over secondary emission");
    }

    string gas = _hasGasEmission ? " and gas emission" : "";
    if (_hasStochasticDustEmission)
        log->info("  Including stochastic dust emission" + gas);
    else if (_hasDustEmission)
        log->info("  Including equilibrium dust emission" + gas);
    else if (_hasGasEmission)
        log->info("  Including gas emission");

    if (_hasSemiDynamicState) log->info("  With semi-dynamic state");

    if (_hasPolarization) log->info("  Including support for polarization");
    if (_hasMovingMedia) log->info("  Including support for kinematics");

    // --- log cosmology ---

    if (_redshift)
    {
        log->info("  Redshift: " + StringUtils::toString(_redshift, 'g', 6));
        double dL = _luminosityDistance / Constants::pc() / 1e6;
        log->info("  Luminosity distance: " + StringUtils::toString(dL, 'g', 6) + " Mpc");
    }

    // --- log model symmetries ---

    // if there are no media, simply log the source model symmetry
    // if there are media, compare the model symmetry to the grid symmetry
    // (the case where the grid has insufficient dimension causes a fatal error in setupSelfBefore)
    if (!_hasMedium)
    {
        log->info("  Model symmetry: " + std::to_string(_modelDimension) + "D");
    }
    else if (_modelDimension == _gridDimension)
    {
        log->info("  Model and grid symmetry: " + std::to_string(_modelDimension) + "D");
    }
    else
    {
        log->info("  Model symmetry: " + std::to_string(_modelDimension)
                  + "D; Spatial grid symmetry: " + std::to_string(_gridDimension) + "D");
        log->warning("  Selecting a grid with the model symmetry might be more efficient");
    }

    // --- other warnings ---

    // enable forced scattering when we have a radiation field because
    // the photon cycle without forced scattering does not support storing the radiation field
    if (_hasRadiationField && !_forceScattering)
    {
        log->warning("  Enabling forced scattering to allow storing the radiation field");
        _forceScattering = true;
    }

    // disable path length stretching if the wavelength of a photon packet can change during its lifetime
    if ((_hasMovingMedia || _hubbleExpansionRate || _hasLymanAlpha) && _forceScattering && _pathLengthBias > 0.)
    {
        log->warning("  Disabling path length stretching to allow Doppler shifts to be properly sampled");
        _pathLengthBias = 0.;
    }

    // if there is a magnetic field, there usually should be spheroidal particles
    if (_hasMagneticField && !_hasSpheroidalPolarization)
        log->warning("  No media have spheroidal particles that could align with the specified magnetic field");
}

////////////////////////////////////////////////////////////////////

void Configuration::setEmulationMode()
{
    _emulationMode = true;
    _hasDynamicState = false;
    _hasPrimaryIterations = false;
    _hasSecondaryIterations = false;
    _includePrimaryEmission = false;
    _minPrimaryIterations = 0;
    _maxPrimaryIterations = 0;
    _minSecondaryIterations = 0;
    _maxSecondaryIterations = 0;
    _numPrimaryPackets = 0;
    _numPrimaryIterationPackets = 0;
    _numSecondaryPackets = 0;
    _numSecondaryIterationPackets = 0;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // This function extends the specified wavelength range with the range of the specified wavelength grid
    void extendForWavelengthGrid(Range& range, WavelengthGrid* wavelengthGrid)
    {
        // we explicitly call setup() on wavelength grids before accessing them
        // because this function may be called early during simulation setup
        wavelengthGrid->setup();
        range.extend(wavelengthGrid->wavelengthRange());
    }

    // This function extends the specified wavelength range with the ranges requested by simulation items
    // in the specified hierarchy that implement the MaterialWavelengthRangeInterface.
    // The function calls itself recursively.
    void extendForMaterialWavelengthRange(Range& range, Item* item)
    {
        auto interface = dynamic_cast<MaterialWavelengthRangeInterface*>(item);
        if (interface)
        {
            Range requested = interface->wavelengthRange();
            if (requested.min() > 0) range.extend(requested);
        }
        for (auto child : item->children()) extendForMaterialWavelengthRange(range, child);
    }
}

////////////////////////////////////////////////////////////////////

Range Configuration::simulationWavelengthRange() const
{
    // include primary and secondary source ranges
    Range range = _sourceWavelengthRange;
    if (_dustEmissionWLG) extendForWavelengthGrid(range, _dustEmissionWLG);

    // extend this range with a wide margin for kinematics if needed
    if (_hasMovingSources || _hasMovingMedia) range.extendWithRedshift(1. / 3.);

    // include radiation field wavelength grid (because dust properties are pre-calculated on these wavelengths)
    if (_hasRadiationField)
    {
        extendForWavelengthGrid(range, _radiationFieldWLG);
        // the calculation of the Planck-integrated absorption cross sections needs this wavelength range;
        // see the precalculate() function in EquilibriumDustEmissionCalculator and StochasticDustEmissionCalculator
        range.extend(Range(0.09e-6, 2000e-6));
    }

    // include default instrument wavelength grid
    if (_defaultWavelengthGrid) extendForWavelengthGrid(range, _defaultWavelengthGrid);

    // include instrument-specific wavelength grids
    auto is = find<InstrumentSystem>(false);
    for (auto ins : is->instruments())
    {
        if (ins->wavelengthGrid()) extendForWavelengthGrid(range, ins->wavelengthGrid());
    }

    // include wavelength ranges requested by simulation items that implement the MaterialWavelengthRangeInterface
    auto sim = find<MonteCarloSimulation>(false);
    extendForMaterialWavelengthRange(range, sim);

    // extend the final range with a narrow margin for round-offs
    range.extendWithRedshift(1. / 100.);
    return range;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // This function adds the characteristic wavelengths of the specified grid to the specified set of wavelengths
    void addForWavelengthGrid(std::set<double>& wavelengths, WavelengthGrid* wavelengthGrid)
    {
        // we explicitly call setup() on wavelength grids before accessing them
        // because this function may be called early during simulation setup
        wavelengthGrid->setup();
        int n = wavelengthGrid->numBins();
        for (int ell = 0; ell != n; ++ell) wavelengths.insert(wavelengthGrid->wavelength(ell));
    }

    // This function adds to the specified set of wavelengths any wavelengths requested by simulation items
    // in the specified hierarchy that implement the MaterialWavelengthRangeInterface.
    // The function calls itself recursively.
    void addForMaterialWavelengthRange(std::set<double>& wavelengths, Item* item)
    {
        auto interface = dynamic_cast<MaterialWavelengthRangeInterface*>(item);
        if (interface)
        {
            // if the range indicates a single nonzero wavelength, then add that wavelength
            Range range = interface->wavelengthRange();
            if (range.min() > 0 && range.min() == range.max()) wavelengths.insert(range.min());

            // if there is a wavelength grid, add it as well
            auto grid = interface->materialWavelengthGrid();
            if (grid) addForWavelengthGrid(wavelengths, grid);
        }
        for (auto child : item->children()) addForMaterialWavelengthRange(wavelengths, child);
    }
}

////////////////////////////////////////////////////////////////////

vector<double> Configuration::simulationWavelengths() const
{
    std::set<double> wavelengths;

    // include radiation field and dust emission wavelength grids
    if (_hasRadiationField) addForWavelengthGrid(wavelengths, _radiationFieldWLG);
    if (_dustEmissionWLG) addForWavelengthGrid(wavelengths, _dustEmissionWLG);

    // include default instrument wavelength grid
    if (_defaultWavelengthGrid) addForWavelengthGrid(wavelengths, _defaultWavelengthGrid);

    // include instrument-specific wavelength grids
    auto is = find<InstrumentSystem>(false);
    for (auto ins : is->instruments())
    {
        if (ins->wavelengthGrid()) addForWavelengthGrid(wavelengths, ins->wavelengthGrid());
    }

    // include wavelength ranges requested by simulation items that implement the MaterialWavelengthRangeInterface
    auto sim = find<MonteCarloSimulation>(false);
    addForMaterialWavelengthRange(wavelengths, sim);

    return vector<double>(wavelengths.begin(), wavelengths.end());
}

////////////////////////////////////////////////////////////////////

WavelengthGrid* Configuration::wavelengthGrid(WavelengthGrid* localWavelengthGrid) const
{
    auto result = localWavelengthGrid && !_oligochromatic ? localWavelengthGrid : _defaultWavelengthGrid;
    if (!result) throw FATALERROR("Cannot find a wavelength grid for instrument or probe");
    return result;
}

////////////////////////////////////////////////////////////////////
