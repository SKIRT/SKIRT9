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
    configuration options, including those offered by all members of the SimulationMode class
    hierarchy.

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
        selected SimulationMode subclass. If any conflicts are found, the function throws a fatal
        error. */
    void setupSelfBefore() override;

    /** This function logs some aspects of the configuration as information to the user. */
    void setupSelfAfter() override;

    //======== Setters that override the user configuration =======

public:
    /** This function puts the simulation in emulation mode. Specifically, it sets a flag that can
        be queried by other simulation items, it sets the number of photon packets to zero, and if
        iteration over the simulation state is enabled, it forces the number of iterations to one.
        */
    void setEmulationMode();

    //=========== Getters for configuration properties ============

public:
    /** Returns true if the simulation has been put in emulation mode. */
    bool emulationMode() const { return _emulationMode; }

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

    /** Returns the number of photon packets launched per primary emission simulation segment. */
    double numPrimaryPackets() const { return _numPrimaryPackets; }

    /** Returns the number of photon packets launched per secondary emission simulation segment. */
    double numIterationPackets() const { return _numIterationPackets; }

    /** Returns the number of photon packets launched per secondary emission simulation segment. */
    double numSecondaryPackets() const { return _numSecondaryPackets; }

    /** Returns true if there is at least one medium component in the simulation. */
    bool hasMedium() const { return _hasMedium; }

    /** Returns the minimum weight reduction factor before a photon packet is terminated. */
    double minWeightReduction() const { return _minWeightReduction; }

    /** Returns the minimum number of forced scattering events before a photon packet is
        terminated. */
    int minScattEvents() const { return _minScattEvents; }

    /** Returns the fraction of path lengths sampled from a linear rather than an exponential
        distribution. */
    double pathLengthBias() const { return _pathLengthBias; }

    /** Returns the number of random density samples for determining spatial cell mass. */
    int numDensitySamples() const { return _numDensitySamples; }

    /** Returns true if the radiation field must be stored during the photon cycle, and false otherwise. */
    bool hasRadiationField() const { return _hasRadiationField; }

    /** Returns true if a panchromatic radiation field (from which a temperature can be calculated)
        is being stored during the photon cycle, and false otherwise. */
    bool hasPanRadiationField() const { return _hasPanRadiationField; }

    /** Returns true if the radiation field for emission from secondary sources must be stored
        (in a separate data structure), and false otherwise. */
    bool hasSecondaryRadiationField() const { return _hasSecondaryRadiationField; }

    /** Returns true if secondary emission must be calculated for any media type, and false otherwise. */
    bool hasSecondaryEmission() const { return _hasDustEmission; }

    /** Returns true if secondary dust emission must be calculated, and false otherwise. */
    bool hasDustEmission() const { return _hasDustEmission; }

    /** Returns true if secondary dust emission must be calculated by taking stochastically heated
        grains into account, and false otherwise. */
    bool hasStochasticDustEmission() const { return _hasStochasticDustEmission; }

    /** Returns true if the cosmic microwave background (CMB) must be added as a source term for
        dust heating, and false otherwise. */
    bool includeHeatingByCMB() const { return _includeHeatingByCMB; }

    /** Returns true if dust self-absorption must be self-consistently calculated through
        iteration, and false otherwise. */
    bool hasDustSelfAbsorption() const { return _hasDustSelfAbsorption; }

    /** Return true if the opacity must be self-consistently calculated through iteration, and
        false otherwise. */
    bool hasOpacityIteration() const { return _hasOpacityIteration; }

    /** Returns the wavelength grid to be used for storing the radiation field. */
    DisjointWavelengthGrid* radiationFieldWLG() const { return _radiationFieldWLG; }

    /** Returns the wavelength grid to be used for calculating the dust emission spectrum. */
    DisjointWavelengthGrid* dustEmissionWLG() const { return _dustEmissionWLG; }

    /** Returns true if the radiation field must be stored during emission (for probing), and false
        otherwise. */
    bool storeEmissionRadiationField() const { return _storeEmissionRadiationField; }

    /** Returns the cell library mapping to be used for calculating the dust emission spectra. */
    SpatialCellLibrary* cellLibrary() const { return _cellLibrary; }

    /** Returns the fraction of secondary photon packets distributed uniformly across spatial
        cells. */
    double secondarySpatialBias() const { return _secondarySpatialBias; }

    /** Returns the fraction of secondary photon packet wavelengths sampled from a bias
        distribution. */
    double secondaryWavelengthBias() const { return _secondaryWavelengthBias; }

    /** Returns the bias distribution for sampling secondary photon packet wavelengths. */
    WavelengthDistribution* secondaryWavelengthBiasDistribution() const { return _secondaryWavelengthBiasDistribution; }

    /** Returns the minimum number of self-absorption iterations. */
    int minIterations() const { return _minIterations; }

    /** Returns the maximum number of self-absorption iterations. */
    int maxIterations() const { return _maxIterations; }

    /** Returns the self-absorption iteration convergece criterion described as follows:
        convergence is reached when the total absorbed dust luminosity is less than this fraction
        of the total absorbed primary luminosity. */
    double maxFractionOfPrimary() const { return _maxFractionOfPrimary; }

    /** Returns the self-absorption iteration convergece criterion described as follows:
        convergence is reached when the total absorbed dust luminosity has changed by less than
        this fraction compared to the previous iteration. */
    double maxFractionOfPrevious() const { return _maxFractionOfPrevious; }

    /** Returns the symmetry dimension of the input model, including sources and media, if present.
        A value of 1 means spherical symmetry, 2 means axial symmetry and 3 means none of these
        symmetries. */
    int modelDimension() const { return _modelDimension; }

    /** Returns the symmetry dimension of the spatial grid, if present, or 0 if there is no spatial
        grid (which can only happen if the simulation does not include any media). A value of 1
        means spherical symmetry, 2 means axial symmetry and 3 means none of these symmetries. */
    int gridDimension() const { return _gridDimension; }

    /** Returns true if at least one medium component in the simulation can have a nonzero velocity
        for some positions, causing Doppler shifts in the perceived, scattered and emitted
        wavelengths. If the function returns false, none of the media has a velocity. */
    bool hasMovingMedia() const { return _hasMovingMedia; }

    /** Returns true if the material mix for at least one medium component in the simulation may
        vary depending on spatial position. If the function returns false, the material mixes and
        thus the material properties for all media are constant throughout the complete spatial
        domain of the simulation. */
    bool hasVariableMedia() const { return _hasVariableMedia; }

    /** Returns true if all media in the simulation support polarization, and false if none of the
        media do. A mixture of support and no support for polarization is not allowed and will
        cause a fatal error during setup. */
    bool hasPolarization() const { return _hasPolarization; }

    /** Returns true if some of the media in the simulation represent spheroidal (i.e.
        non-spherical) particles and require the corresponding treatment of polarization for
        scattering, absorption and emission, or false otherwise. If this function returns true, the
        hasPolarization() and hasMagneticField() functions return true as well. */
    bool hasSpheroidalPolarization() const { return _hasSpheroidalPolarization; }

    /** Returns true if a medium component in the simulation defines a spatial magnetic field
        distribution that may have nonzero strength for some positions, or false if none of the
        media define a magnetic field. It is not allowed for multiple medium components to define
        a magnetic field (a fatal error is raised during setup when this happens). */
    bool hasMagneticField() const { return _hasMagneticField; }

    //======================== Data Members ========================

private:
    // general
    bool _emulationMode{false};

    // cosmology parameters
    double _redshift{0.};
    double _angularDiameterDistance{0.};
    double _luminosityDistance{0.};

    // primary source wavelengths
    bool _oligochromatic{false};
    Range _sourceWavelengthRange;
    WavelengthGrid* _defaultWavelengthGrid{nullptr};
    WavelengthDistribution* _oligoWavelengthBiasDistribution{nullptr};

    // launch
    double _numPrimaryPackets{0.};
    double _numIterationPackets{0.};
    double _numSecondaryPackets{0.};

    // extinction
    bool _hasMedium{false};
    double _minWeightReduction{1e4};
    int _minScattEvents{0};
    double _pathLengthBias{0.5};
    int _numDensitySamples{100};

    // radiation field
    bool _hasRadiationField{false};
    bool _hasPanRadiationField{false};
    bool _hasSecondaryRadiationField{false};
    DisjointWavelengthGrid* _radiationFieldWLG{nullptr};

    // emission
    bool _hasDustEmission{false};
    bool _hasStochasticDustEmission{false};
    bool _includeHeatingByCMB{false};
    bool _hasDustSelfAbsorption{false};
    DisjointWavelengthGrid* _dustEmissionWLG{nullptr};
    SpatialCellLibrary* _cellLibrary{nullptr};
    bool _storeEmissionRadiationField{false};
    double _secondarySpatialBias{0.5};
    double _secondaryWavelengthBias{0.5};
    WavelengthDistribution* _secondaryWavelengthBiasDistribution{nullptr};
    int _minIterations{1};
    int _maxIterations{10};
    double _maxFractionOfPrimary{0.01};
    double _maxFractionOfPrevious{0.03};

    // opacity iteration
    bool _hasOpacityIteration{false};

    // properties derived from the configuration at large
    int _modelDimension{0};
    int _gridDimension{0};
    bool _hasMovingSources{false};
    bool _hasMovingMedia{false};
    bool _hasVariableMedia{false};
    bool _hasPolarization{false};
    bool _hasSpheroidalPolarization{false};
    bool _hasMagneticField{false};
};

////////////////////////////////////////////////////////////////////

#endif
