/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

#include "Array.hpp"
#include "Range.hpp"
#include "SimulationItem.hpp"
class DisjointWavelengthGrid;
class SpatialCellLibrary;
class WavelengthDistribution;
class WavelengthGrid;

////////////////////////////////////////////////////////////////////

/** Configuration is a helper class that serves as a central clearing house for overall simulation
    configuration options including the simulation mode.

    Each MonteCarloSimulation holds a single Configuration object. During setup, it retrieves many
    properties and options from the simulation hierarchy, verifying consistency of the
    configuration and flagging any conflicts while doing so. Once this process has completed, the
    Configuration object offers getters for these retrieved properties to any of the other
    simulation items in the hierarchy. The setup() function of the Configuration object is invoked
    at the very early stages of simulation setup, so that it is safe for other simulation items to
    retrieve information from the Configuration object during setup.

    The Configuration class is based on SimulationItem so that it can be part of a simulation item
    hierarchy, however it is not discoverable because it is not intended to be selected or
    configured by the user. */
class Configuration : public SimulationItem
{
    //============= Construction - Setup - Destruction =============

public:
    /** This constructor creates a Configuration object that is hooked up as a child to the
        specified parent in the simulation hierarchy, so that it will automatically be deleted. The
        setup() function is \em not called by this constructor. */
    explicit Configuration(SimulationItem* parent);

protected:
    /** This function retrieves properties and options from the simulation hierarchy and stores the
        resulting values internally so that they can be returned by any of the getters with minimal
        overhead. During this process, the function also verifies the consistency of the simulation
        configuration, for example checking the configuration against the requirements of the
        selected simulation mode. If any conflicts are found, the function throws a fatal error. */
    void setupSelfBefore() override;

    /** This function logs some aspects of the configuration as information to the user. */
    void setupSelfAfter() override;

    //======== Setters that override the user configuration =======

public:
    /** This function puts the simulation in emulation mode. Specifically, it sets a flag that can
        be queried by other simulation items, it sets the number of photon packets to zero, and
        disables iteration over primary and/or secondary emisson. */
    void setEmulationMode();

    //=========== Getters for configuration properties ============

public:
    // ----> emulation mode

    /** Returns true if the simulation has been put in emulation mode. */
    bool emulationMode() const { return _emulationMode; }

    // ----> symmetry

    /** Returns the symmetry dimension of the input model, including sources and media, if present.
        A value of 1 means spherical symmetry, 2 means axial symmetry and 3 means none of these
        symmetries. */
    int modelDimension() const { return _modelDimension; }

    /** Returns the symmetry dimension of the spatial grid, if present, or 0 if there is no spatial
        grid (which can only happen if the simulation does not include any media). A value of 1
        means spherical symmetry, 2 means axial symmetry and 3 means none of these symmetries. */
    int gridDimension() const { return _gridDimension; }

    // ----> cosmology

    /** Returns the redshift at which the model resides, or zero if the model resides in the Local
        Universe. */
    double redshift() const { return _redshift; }

    /** Returns the angular-diameter distance corresponding to the redshift at which the model
        resides, or zero if the model resides in the Local Universe. Refer to the Cosmology class
        description for more information. */
    double angularDiameterDistance() const { return _angularDiameterDistance; }

    /** Returns the luminosity distance corresponding to the redshift at which the model resides,
        or zero if the model resides in the Local Universe. Refer to the Cosmology class
        description for more information. */
    double luminosityDistance() const { return _luminosityDistance; }

    // ----> wavelengths

    /** Returns true if the wavelength regime of the simulation is oligochromatic. */
    bool oligochromatic() const { return _oligochromatic; }

    /** Returns the total wavelength range of the primary sources in the simulation. For
        panchromatic simulations, this range is configured by the user in the source system. For
        oligochromatic simulations, the range includes the discrete source wavelengths used in the
        simulation, which are also user-configured in the source system. */
    Range sourceWavelengthRange() const { return _sourceWavelengthRange; }

    /** Returns a wavelength range that covers all wavelengths possibly used in the simulation for
        photon transport or for otherwise probing material properties (e.g. optical depth). This
        range includes the primary and secondary source wavelength ranges extended on both sides to
        accommodate a redshift or blueshift caused by kinematics corresponding to \f$v/c=1/3\f$. It
        also includes the range of the instrument wavelength grids and the wavelengths used for
        material normalization and material property probes. */
    Range simulationWavelengthRange() const;

    /** Returns a list of wavelengths that are explicitly or indirectly mentioned by the simulation
        configuration. This includes the characteristic wavelengths of all configured wavelength
        grids (for instruments, probes, radiation field or dust emission) and specific wavelengths
        used for normalization or probing. */
    vector<double> simulationWavelengths() const;

    /** Returns the wavelength grid to be used for an instrument or probe, given the wavelength
        grid configured locally for the calling instrument or probe (which may the null pointer to
        indicate that no local grid was configured). For oligochromatic simulations, the function
        always returns a wavelength grid with disjoint bins centered around the discrete source
        wavelengths used in the simulation. For panchromatic simulations, the function returns the
        provided local wavelength grid if it is non-null, and otherwise it returns the default
        instrument wavelength grid obtained from the instrument system. If both the provided local
        wavelength grid and the default instrument wavelength grid are the null pointer, the
        function throws a fatal error. */
    WavelengthGrid* wavelengthGrid(WavelengthGrid* localWavelengthGrid) const;

    /** For oligochromatic simulations, this function returns the wavelength bias distribution to
        be used by all primary sources. For panchromatic simulations, the function returns the null
        pointer. */
    WavelengthDistribution* oligoWavelengthBiasDistribution() { return _oligoWavelengthBiasDistribution; }

    // ----> probes

    /** Returns true if one of the Snapshot::getEntities() functions may be called for any of the
        snapshots associated with the imported sources and media in the simulation, implying that
        the snapshot must prebuild the required search data structures. In the current
        implementation, this happens only if the simulation includes one or more input model
        probes, i.e. instances of an InputModelProbe subclass. */
    bool snapshotsNeedGetEntities() const { return _snapshotsNeedGetEntities; }

    // ----> media

    /** Returns true if there is at least one medium component in the simulation, and false
        otherwise. */
    bool hasMedium() const { return _hasMedium; }

    /** Returns true if the Medium::generatePosition() function may be called for the media in the
        simulation. In the current implementation, this happens only if the simulation uses a
        VoronoiMeshSpatialGrid instance to discretize the spatial domain. If there are no media or
        the Medium::generatePosition() will never be called during this simulation, this function
        returns false. */
    bool mediaNeedGeneratePosition() const { return _mediaNeedGeneratePosition; }

    /** Returns true if one or more medium components in the simulation may have a nonzero velocity
        for some positions. If the function returns false, none of the media has a velocity. */
    bool hasMovingMedia() const { return _hasMovingMedia; }

    /** Returns true if one of the medium components in the simulation defines a spatial magnetic
        field distribution that may have nonzero strength for some positions, or false if none of
        the media define a magnetic field. It is not allowed for multiple medium components to
        define a magnetic field (a fatal error is raised during setup when this happens). */
    bool hasMagneticField() const { return _magneticFieldMediumIndex >= 0; }

    /** Returns the index of the medium component defining the magnetic field, or a negative value
        if none of the media define a magnetic field. */
    int magneticFieldMediumIndex() const { return _magneticFieldMediumIndex; }

    /** Returns true if the material mix for at least one medium component in the simulation may
        vary depending on spatial position. If the function returns false, the material mixes and
        thus the material properties for all media are constant throughout the complete spatial
        domain of the simulation. */
    bool hasVariableMedia() const { return _hasVariableMedia; }

    /** Returns true if the perceived photon packet wavelength equals the intrinsic photon packet
        wavelength for all spatial cells along the path of the packet. The following conditions
        cause this function to return false: Hubble expansion is enabled or some media may have a
        non-zero velocity in some cells. */
    bool hasConstantPerceivedWavelength() const { return _hasConstantPerceivedWavelength; }

    /** Returns true if the simulation has a exactly one medium component and the absorption and
        scattering cross sections for a photon packet traversing that medium component are
        spatially constant, so that the opacity in each crossed cell can be calculated by
        multiplying this constant cross section by the number density in the cell. Otherwise the
        function returns false.

        The following conditions cause this function to return false: Hubble expansion is enabled,
        there is more than one medium component, the medium may have a non-zero velocity in some
        cells, the medium has a variable material mix; the cross sections for some material mixes
        depend on extra medium state variables such as temperature or fragment weight factors. */
    bool hasSingleConstantSectionMedium() const { return _hasSingleConstantSectionMedium; }

    /** Returns true if the simulation has two or more medium components and the absorption and
        scattering cross sections for a photon packet traversing those medium components are
        spatially constant, so that the opacity in each crossed cell can be calculated by
        multiplying these constant cross sections by the corresponding number densities in the
        cell. Otherwise the function returns false.

        The following conditions cause this function to return false: Hubble expansion is enabled,
        some media may have a non-zero velocity in some cells, so that the perceived wavelength
        changes between cells; some media have a variable material mix; the cross sections for some
        material mixes depend on extra medium state variables such as temperature or fragment
        weight factors. */
    bool hasMultipleConstantSectionMedia() const { return _hasMultipleConstantSectionMedia; }

    /** Returns true if a scattering interaction for one or more media may emulate secondary
        emission, and false otherwise. */
    bool scatteringEmulatesSecondaryEmission() const { return _scatteringEmulatesSecondaryEmission; }

    /** Returns true if a scattering interaction for one or more media may adjust the wavelength of
        the interacting photon packet or may emulate secondary emission, and false otherwise. */
    bool needIndividualPeelOff() const { return _needIndividualPeelOff; }

    /** Returns true if all media in the simulation support polarization, and false if none of the
        media do. A mixture of support and no support for polarization is not allowed and will
        cause a fatal error during setup. */
    bool hasPolarization() const { return _hasPolarization; }

    /** Returns true if some of the media in the simulation represent spheroidal (i.e.
        non-spherical) particles and require the corresponding treatment of polarization for
        scattering, absorption and emission, or false otherwise. If this function returns true, the
        hasPolarization() and hasMagneticField() functions return true as well. */
    bool hasSpheroidalPolarization() const { return _hasSpheroidalPolarization; }

    // ----> media sampling

    /** Returns the number of random spatial samples for determining density (or mass). */
    int numDensitySamples() const { return _numDensitySamples; }

    /** Returns the number of random spatial samples for determining other properties. */
    int numPropertySamples() const { return _numPropertySamples; }

    // ----> phases, iterations, number of packets

    /** Returns true if secondary emission must be calculated for any media type, and false
        otherwise. */
    bool hasSecondaryEmission() const { return _hasSecondaryEmission; }

    /** Returns true if the simulation iterates over primary emission, and false otherwise. */
    bool hasPrimaryIterations() const { return _hasPrimaryIterations; }

    /** Returns true if the simulation iterates over secondary emission, and false otherwise. */
    bool hasSecondaryIterations() const { return _hasSecondaryIterations; }

    /** Returns true if the simulation iterates over both primary and secondary emission and the
        iterations over secondary emission should include primary emission, and false otherwise. */
    bool hasMergedIterations() const { return _hasMergedIterations; }

    /** Returns the minimum number of iterations in the primary emission phase. */
    int minPrimaryIterations() const { return _minPrimaryIterations; }

    /** Returns the maximum number of iterations in the primary emission phase. */
    int maxPrimaryIterations() const { return _maxPrimaryIterations; }

    /** Returns the minimum number of iterations in the secondary emission phase. */
    int minSecondaryIterations() const { return _minSecondaryIterations; }

    /** Returns the maximum number of iterations in the secondary emission phase. */
    int maxSecondaryIterations() const { return _maxSecondaryIterations; }

    /** Returns the number of photon packets launched per regular primary emission simulation
        segment. */
    double numPrimaryPackets() const { return _numPrimaryPackets; }

    /** Returns the number of photon packets launched per iteration segment during primary
        emission. */
    double numPrimaryIterationPackets() const { return _numPrimaryIterationPackets; }

    /** Returns the number of photon packets launched per regular secondary emission simulation
        segment. */
    double numSecondaryPackets() const { return _numSecondaryPackets; }

    /** Returns the number of photon packets launched per iteration segment during secondary
        emission. */
    double numSecondaryIterationPackets() const { return _numSecondaryIterationPackets; }

    /** Returns the dust self-absorption iteration convergence criterion described as follows:
        convergence is reached when the total absorbed dust luminosity is less than this fraction
        of the total absorbed primary luminosity. */
    double maxFractionOfPrimary() const { return _maxFractionOfPrimary; }

    /** Returns the dust self-absorption iteration convergence criterion described as follows:
        convergence is reached when the total absorbed dust luminosity has changed by less than
        this fraction compared to the previous iteration. */
    double maxFractionOfPrevious() const { return _maxFractionOfPrevious; }

    // ----> dynamic medium state

    /** Returns true if the simulation has primary or merged iterations and includes one or more
        dynamic medium state recipes (instances of a DynamicStateRecipe subclass), and false
        otherwise. In the current implementation, these recipes by definition perform primary
        dynamic medium state (PDMS) updates. */
    bool hasDynamicStateRecipes() const { return _hasDynamicStateRecipes; }

    /** Returns true if the simulation has primary or merged iterations and includes one or more
        media with an associated MaterialMix that performs primary dynamic medium state (PDMS)
        updates, and false otherwise. */
    bool hasPrimaryDynamicStateMedia() const { return _hasPrimaryDynamicStateMedia; }

    /** Returns true if the simulation has secondary emission and includes one or more media with
        an associated MaterialMix that performs secondary dynamic medium state (SDMS) updates, and
        false otherwise. */
    bool hasSecondaryDynamicStateMedia() const { return _hasSecondaryDynamicStateMedia; }

    /** Returns true if the simulation has primary or merged iterations and includes one or more
        recipes or media that perform primary dynamic medium state (PDMS) updates, and false
        otherwise. */
    bool hasPrimaryDynamicState() const { return _hasPrimaryDynamicState; }

    /** Returns true if the simulation has primary or merged iterations and includes one or more
        recipes or media that perform secondary dynamic medium state (SDMS) updates, and false
        otherwise. */
    bool hasSecondaryDynamicState() const { return _hasSecondaryDynamicState; }

    /** Returns true if the simulation has primary or merged iterations and includes one or more
        recipes or media that perform primary or secondary dynamic medium state (DMS) updates, and
        false otherwise. */
    bool hasDynamicState() const { return _hasPrimaryDynamicState || _hasSecondaryDynamicState; }

    // ----> photon cycle

    /** Returns true if the extinction cross section (the sum of the absorption and scattering
        cross section) for one or more material mixes in the simulation can be negative, and false
        if not. */
    bool hasNegativeExtinction() const { return _hasNegativeExtinction; }

    /** Returns true if explicit absorption should be used during the photon cycle, false if not.
        */
    bool explicitAbsorption() const { return _explicitAbsorption; }

    /** Returns true if forced scattering should be used during the photon cycle, false if not. */
    bool forceScattering() const { return _forceScattering; }

    /** Returns the minimum weight reduction factor before a photon packet is terminated. */
    double minWeightReduction() const { return _minWeightReduction; }

    /** Returns the minimum number of forced scattering events before a photon packet is
        terminated. */
    int minScattEvents() const { return _minScattEvents; }

    /** Returns the fraction of path lengths sampled from a linear rather than an exponential
        distribution. */
    double pathLengthBias() const { return _pathLengthBias; }

    /** This enumeration lists the supported Lyman-alpha acceleration schemes. */
    enum class LyaAccelerationScheme { None, Constant, Variable };

    /** Returns the enumeration value determining the acceleration scheme to be used for
        Lyman-alpha line scattering. The value is relevant only if Lyman-alpha line treatment is
        enabled in the simulation. */
    LyaAccelerationScheme lyaAccelerationScheme() const { return _lyaAccelerationScheme; }

    /** Returns the strength of the Lyman-alpha acceleration scheme to be applied. The value is
        relevant only if Lyman-alpha line treatment is enabled in the simulation and
        lyaAccelerationScheme() returns \c Constant or \c Variable. */
    double lyaAccelerationStrength() const { return _lyaAccelerationStrength; }

    /** If inclusion of the Hubble flow is enabled, this function returns the relative expansion
        rate of the universe in which the model resides. If inclusion of the Hubble flow is
        disabled, or if the simulation does not include Lyman-alpha treatment, this function
        returns zero. */
    double hubbleExpansionRate() const { return _hubbleExpansionRate; }

    // ----> radiation field

    /** Returns true if the radiation field must be stored during the photon cycle, and false
        otherwise. */
    bool hasRadiationField() const { return _hasRadiationField; }

    /** Returns true if a panchromatic radiation field (from which a temperature can be calculated)
        is being stored during the photon cycle, and false otherwise. */
    bool hasPanRadiationField() const { return _hasPanRadiationField; }

    /** Returns true if the radiation field for emission from secondary sources must be stored in a
        separate data structure, and false otherwise. */
    bool hasSecondaryRadiationField() const { return _hasSecondaryRadiationField; }

    /** Returns the wavelength grid to be used for storing the radiation field, or the null pointer
        if hasRadiationField() returns false. */
    DisjointWavelengthGrid* radiationFieldWLG() const { return _radiationFieldWLG; }

    // ----> secondary emission

    /** Returns true if the radiation field must be stored during emission (for probing), and false
        otherwise. */
    bool storeEmissionRadiationField() const { return _storeEmissionRadiationField; }

    /** Returns the fraction of secondary photon packets distributed uniformly across spatial
        cells. */
    double secondarySpatialBias() const { return _secondarySpatialBias; }

    /** Returns the fraction of secondary photon packets distributed uniformly across secondary
        sources. */
    double secondarySourceBias() const { return _secondarySourceBias; }

    // ----> dust emission

    /** Returns true if thermal dust emission must be calculated, and false otherwise. */
    bool hasDustEmission() const { return _hasDustEmission; }

    /** Returns true if thermal dust emission must be calculated by taking stochastically heated
        grains into account, and false otherwise. */
    bool hasStochasticDustEmission() const { return _hasStochasticDustEmission; }

    /** Returns true if the cosmic microwave background (CMB) must be added as a source term for
        dust heating, and false otherwise. */
    bool includeHeatingByCMB() const { return _includeHeatingByCMB; }

    /** Returns the wavelength grid to be used for calculating the dust emission spectra. */
    DisjointWavelengthGrid* dustEmissionWLG() const { return _dustEmissionWLG; }

    /** Returns the cell library mapping to be used for calculating the dust emission spectra. */
    SpatialCellLibrary* cellLibrary() const { return _cellLibrary; }

    /** Returns the bias weight for dust emission sources. */
    double dustEmissionSourceWeight() const { return _dustEmissionSourceWeight; }

    /** Returns the fraction of dust emission photon packet wavelengths sampled from a bias
        distribution. */
    double dustEmissionWavelengthBias() const { return _dustEmissionWavelengthBias; }

    /** Returns the bias distribution for sampling dust emission photon packet wavelengths. */
    WavelengthDistribution* dustEmissionWavelengthBiasDistribution() const
    {
        return _dustEmissionWavelengthBiasDistribution;
    }

    // ----> gas emission

    /** Returns true if gas emission must be calculated, and false otherwise. */
    bool hasGasEmission() const { return _hasGasEmission; }

    //======================== Data Members ========================

private:
    // emulation mode
    bool _emulationMode{false};

    // symmetry
    int _modelDimension{0};
    int _gridDimension{0};

    // cosmology
    double _redshift{0.};
    double _angularDiameterDistance{0.};
    double _luminosityDistance{0.};

    // wavelengths
    bool _oligochromatic{false};
    Range _sourceWavelengthRange;
    WavelengthGrid* _defaultWavelengthGrid{nullptr};
    WavelengthDistribution* _oligoWavelengthBiasDistribution{nullptr};

    // probes
    bool _snapshotsNeedGetEntities{false};

    // media
    bool _hasMedium{false};
    bool _mediaNeedGeneratePosition{false};
    bool _hasMovingSources{false};
    bool _hasMovingMedia{false};
    int _magneticFieldMediumIndex{-1};
    bool _hasVariableMedia{false};
    bool _hasConstantPerceivedWavelength{false};
    bool _hasSingleConstantSectionMedium{false};
    bool _hasMultipleConstantSectionMedia{false};
    bool _hasScatteringDispersion{false};
    bool _scatteringEmulatesSecondaryEmission{false};
    bool _needIndividualPeelOff{false};
    bool _hasPolarization{false};
    bool _hasSpheroidalPolarization{false};

    // media sampling
    int _numDensitySamples{100};
    int _numPropertySamples{1};

    // phases, iterations, number of packets
    bool _hasSecondaryEmission{false};
    bool _hasPrimaryIterations{false};
    bool _hasSecondaryIterations{false};
    bool _hasMergedIterations{false};
    int _minPrimaryIterations{1};
    int _maxPrimaryIterations{10};
    int _minSecondaryIterations{1};
    int _maxSecondaryIterations{10};
    double _numPrimaryPackets{0.};
    double _numPrimaryIterationPackets{0.};
    double _numSecondaryPackets{0.};
    double _numSecondaryIterationPackets{0.};
    double _maxFractionOfPrimary{0.01};
    double _maxFractionOfPrevious{0.03};

    // dynamic medium state
    bool _hasDynamicStateRecipes{false};
    bool _hasPrimaryDynamicStateMedia{false};
    bool _hasSecondaryDynamicStateMedia{false};
    bool _hasPrimaryDynamicState{false};
    bool _hasSecondaryDynamicState{false};

    // photon cycle
    bool _hasNegativeExtinction{false};
    bool _explicitAbsorption{false};
    bool _forceScattering{true};
    double _minWeightReduction{1e4};
    int _minScattEvents{0};
    double _pathLengthBias{0.5};
    bool _hasLymanAlpha{false};
    LyaAccelerationScheme _lyaAccelerationScheme{LyaAccelerationScheme::Variable};
    double _lyaAccelerationStrength{1.};
    double _hubbleExpansionRate{0.};

    // radiation field
    bool _hasRadiationField{false};
    bool _hasPanRadiationField{false};
    bool _hasSecondaryRadiationField{false};
    DisjointWavelengthGrid* _radiationFieldWLG{nullptr};

    // secondary emission
    bool _storeEmissionRadiationField{false};
    double _secondarySpatialBias{0.5};
    double _secondarySourceBias{0.5};

    // dust emission
    bool _hasDustEmission{false};
    bool _hasStochasticDustEmission{false};
    bool _includeHeatingByCMB{false};
    DisjointWavelengthGrid* _dustEmissionWLG{nullptr};
    SpatialCellLibrary* _cellLibrary{nullptr};
    double _dustEmissionSourceWeight{1.};
    double _dustEmissionWavelengthBias{0.5};
    WavelengthDistribution* _dustEmissionWavelengthBiasDistribution{nullptr};

    // gas emission
    bool _hasGasEmission{false};
};

////////////////////////////////////////////////////////////////////

#endif
