/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FRAGMENTDUSTMIXDECORATOR_HPP
#define FRAGMENTDUSTMIXDECORATOR_HPP

#include "MultiGrainDustMix.hpp"
class DustMixFragment;

////////////////////////////////////////////////////////////////////

/** A FragmentDustMixDecorator instance aggregates a fixed number of predefined dust populations,
    called \em fragments, and provides a way to manage the relative weight of each of those
    fragments during the simulation. The aggregated fragments are defined by specifying a
    MultiGrainDustMix instance as part of the user configuration. Depending on a user configuration
    flag, the specified multi-grain dust mix is fragmented into its underlying dust populations or
    into the individual grain size bins for those populations. A fragment dust mix decorator thus
    essentially splits a multi-grain dust mix into its underlying fragments. The dust mix being
    fragmented can be a turn-key multi-grain mix (e.g., ZubkoDustMix) or a fully user-configured
    dust mix (ConfigurableDustMix), offering nearly limitless possibilities.

    The specific medium state for a fragmented dust mix consists of a weight factor for each
    fragment. To construct the aggregated dust mix represented by the fragmented dust mix for a
    given spatial cell, the dust mass represented by each fragment in the original dust mix is
    multiplied by its corresponding weight factor in the specific state. In other words, a weight
    factor of one means that the fragment's mass remains unchanged from the original dust mix.

    Initialization of the specific medium state proceeds as follows. If the fragmented dust mix is
    configured as part of a geometric medium component, all fragment weight factors are set to a
    value of one, thus preserving the relative weigths of the fragments in the original dust mix.
    If the fragmented dust mix is configured as part of an imported medium component, the fragment
    weight factors are imported from the snapshot.

    FragmentDustMixDecorator inherits directly from MaterialMix rather than from DustMix or
    MultiGrainDustMix so that it cannot be nested inside another fragment dust mix decorator or
    inside a SelectDustMixFamily instance. */
class FragmentDustMixDecorator : public MaterialMix, public MultiGrainPopulationInterface
{
    ITEM_CONCRETE(FragmentDustMixDecorator, MaterialMix,
                  "a dust mix decorator that manages separate densities for fragments of another dust mix")
        ATTRIBUTE_TYPE_INSERT(FragmentDustMixDecorator, "CustomMediumState")
        ATTRIBUTE_TYPE_DISPLAYED_IF(FragmentDustMixDecorator, "Level2")

        PROPERTY_ITEM(dustMix, MultiGrainDustMix, "the dust mix to be fragmented")
        ATTRIBUTE_DEFAULT_VALUE(dustMix, "ThemisDustMix")

        PROPERTY_BOOL(fragmentSizeBins, "fragment each dust population into its grain size bins")
        ATTRIBUTE_DEFAULT_VALUE(fragmentSizeBins, "false")

        PROPERTY_BOOL(hasDynamicDensities, "allow the fragment densities to be adjusted dynamically")
        ATTRIBUTE_RELEVANT_IF(hasDynamicDensities, "DynamicState")
        ATTRIBUTE_DEFAULT_VALUE(hasDynamicDensities, "DynamicState:true;false")
        ATTRIBUTE_DISPLAYED_IF(hasDynamicDensities, "Level3")

        PROPERTY_DOUBLE(initialDensityFraction, "the initial value of the dynamic density fraction")
        ATTRIBUTE_MIN_VALUE(initialDensityFraction, "[0")
        ATTRIBUTE_MAX_VALUE(initialDensityFraction, "1]")
        ATTRIBUTE_DEFAULT_VALUE(initialDensityFraction, "0")
        ATTRIBUTE_RELEVANT_IF(initialDensityFraction, "DynamicState&hasDynamicDensities")
        ATTRIBUTE_DISPLAYED_IF(initialDensityFraction, "Level3")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function creates and remembers a dust mix fragment for each grain population in the
        multi-grain dust mix being fragmented or, if requested by the user, for each size bin in
        each of those grain populations. */
    void setupSelfAfter() override;

    //======== Material type =======

public:
    /** This function returns the fundamental material type represented by this material mix. For
        the fragment dust mix decorator, it returns MaterialType::Dust. */
    MaterialType materialType() const override;

    //======================== Capabilities =======================

public:
    /** This function returns a flag indicating whether the material mix supports polarization
        during scattering events or not. For the fragment dust mix decorator, the function returns
        the value returned by the corresponding function for the dust mix being fragmented. */
    bool hasPolarizedScattering() const override;

    /** This function returns a flag indicating whether this material mix supports stochastic
        heating of dust grains. For the fragment dust mix decorator, the function returns true
        because all dust mixes being fragmented support stochastic heating by definition. */
    bool hasStochasticDustEmission() const override;

    /** This function returns true, indicating that the cross sections returned by the fragment
        dust mix decorator depend on the values of specific state variables other than the number
        density, in this case the fragment weight factors. */
    bool hasExtraSpecificState() const override;

    /** This function returns true for this class because all dust mixes support secondary
        continuum emission. */
    bool hasContinuumEmission() const override;

    //======== Medium state setup =======

public:
    /** This function returns the number and type of import parameters required by this fragmented
        dust mix as a list of SnapshotParameter objects. Specifically, the function returns a
        parameter description for each fragment corresponding to the weight factor to be assigned
        to that fragment. The parameters imported for each spatial cell end up in the state for the
        medium component associated with this fragmented dust mix. For more information, see the
        specificStateVariableInfo() and initializeSpecificState() functions. */
    vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns a list of StateVariable objects describing the specific state
        variables used by this fragmented dust mix. Specifically, in addition to a description for
        the overall number density, the function returns a custom state variable description for
        each fragment corresponding to the weight factor for that fragment. For more information,
        see the parameterInfo() and initializeSpecificState() functions. */
    vector<StateVariable> specificStateVariableInfo() const override;

    /** This function initializes the specific state variables requested by this fragmented dust
        mix through the specificStateVariableInfo() function except for the number density. If the
        \em params array is nonempty, the function assumes that the array contains the appropriate
        number of imported parameter values as requested through the parameterInfo() function, i.e.
        a weight factor per fragment. It then copies these weights into the corresponding specific
        state variables. If the \em params array is empty, the function sets all the weights in the
        specific state to a value of one. The \em metallicity and \em temperature arguments are
        ignored. */
    void initializeSpecificState(MaterialState* state, double metallicity, double temperature,
                                 const Array& params) const override;

    //======== Low-level material properties =======

public:
    /** This function returns the dust mass \f$\mu\f$ per hydrogen atom for the original dust mix
        being fragmented, without taking into account any additional fragment weights. */
    double mass() const override;

    /** This function returns the absorption cross section per entity
        \f$\varsigma^{\text{abs}}_{\lambda}\f$ at wavelength \f$\lambda\f$ for the original dust
        mix being fragmented, without taking into account any additional fragment weights. */
    double sectionAbs(double lambda) const override;

    /** This function returns the scattering cross section per entity
        \f$\varsigma^{\text{sca}}_{\lambda}\f$ at wavelength \f$\lambda\f$ for the original dust
        mix being fragmented, without taking into account any additional fragment weights. */
    double sectionSca(double lambda) const override;

    /** This function returns the total extinction cross section per entity
        \f$\varsigma^{\text{ext}}_{\lambda} = \varsigma^{\text{abs}}_{\lambda} +
        \varsigma^{\text{sca}}_{\lambda}\f$ at wavelength \f$\lambda\f$ for the original dust mix
        being fragmented, without taking into account any additional fragment weights. */
    double sectionExt(double lambda) const override;

    /** This function returns the scattering asymmetry parameter \f$g_\lambda =
        \left<\cos\theta\right>\f$ at wavelength \f$\lambda\f$ for the original dust mix being
        fragmented, without taking into account any additional fragment weights. */
    double asymmpar(double lambda) const override;

    //======== High-level photon life cycle =======

    /** This function returns the absorption opacity \f$k^\text{abs}=n\varsigma^\text{abs}\f$ for
        the given wavelength and material state, taking into account the fragment weights in the
        material state as described in the class header. The photon properties are used only if
        they are used by the dust mix being fragmented. */
    double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the scattering opacity \f$k^\text{sca}=n\varsigma^\text{sca}\f$ for
        the given wavelength and material state, taking into account the fragment weights in the
        material state as described in the class header. The photon properties are used only if
        they are used by the dust mix being fragmented. */
    double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the extinction opacity \f$k^\text{ext}=k^\text{abs}+k^\text{sca}\f$
        for the given wavelength and material state, taking into account the fragment weights in
        the material state as described in the class header. The photon properties are used only if
        they are used by the dust mix being fragmented. */
    double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function calculates the contribution of the medium component associated with this
        material mix to the peel-off photon luminosity, polarization state, and wavelength shift
        for the given wavelength, geometry, material state, and photon properties. The relative
        weight of each fragment's contribution is adjusted by the relative scattering opacity of
        the fragment. */
    void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs, Direction bfky,
                           const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function performs a scattering event on the specified photon packet in the spatial
        cell and medium component represented by the specified material state. Specifically, it
        randomly selects one of the fragments in the fragmented dust mix, weighted by their
        relative scattering opacity, and then causes this fragment to perform the scattering event.
        */
    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======== Secondary emission =======

public:
    /** This function returns the wavelength grid on which dust emission is discretized, i.e. the
        wavelength grid returned by the Configuration::dustEmissionWLG() function. */
    DisjointWavelengthGrid* emissionWavelengthGrid() const override;

    /** This function returns the emissivity spectrum \f$\varepsilon_{\ell'}\f$ (radiated power per
        unit of solid angle and per hydrogen atom) of the dust mix when it would be embedded in a
        given radiation field. Because this function does not have access to the material state, it
        simply returns the value returned by the corresponding function for the dust mix being
        fragmented. In other words, it behaves as if the fragment weights are all equal to one.

        The input and output arrays are discretized on the wavelength grids returned by the
        Configuration::radiationFieldWLG() and Configuration::dustEmissionWLG() functions,
        repectively. */
    Array emissivity(const Array& Jv) const override;

    /** This function returns the emission spectrum (radiated power per unit of solid angle) in the
        spatial cell and medium component represented by the specified material state and the
        receiving material mix when it would be embedded in the specified radiation field. The
        returned spectrum takes into account the weigths of the various fragments in the mix in
        addition to the number density of the material in the specified cell.

        The input and output arrays are discretized on the wavelength grids returned by the
        Configuration::radiationFieldWLG() and Configuration::dustEmissionWLG() functions,
        repectively. */
    Array emissionSpectrum(const MaterialState* state, const Array& Jv) const override;

    /** This function returns an indicative temperature of the fragmented dust mix when it would be
        embedded in a given radiation field, averaging over the fragments using the fragment
        weights in the material state as described in the class header. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //=========== Exposing fragments (MultiGrainPopulationInterface) ============

public:
    /** This function returns the number of fragments (with indices \f$f\f$) in this fragmented
        dust mix. If the \em fragmentSizeBins flag is configured to false, this corresponds to the
        number of populations in the dust mix being fragmented. If the \em fragmentSizeBins flag is
        configured to true, this corresponds to the grand total of size bins for all populations in
        the dust mix being fragmented. */
    int numPopulations() const override;

    /** This function returns a brief human-readable identifier for the type of grain material
        represented by the fragment with index \f$f\f$. */
    string populationGrainType(int f) const override;

    /** This function returns the bulk mass density \f$\rho_\text{bulk}\f$ of the grain material
        represented by the fragment with index \f$f\f$. */
    double populationBulkDensity(int f) const override;

    /** This function returns the minimum and maximum grain sizes \f$a_{\text{min},f},
        a_{\text{max},c}\f$ for the fragment with index \f$f\f$. */
    Range populationSizeRange(int f) const override;

    /** This function returns the grain size distribution object for the fragment with index
        \f$f\f$. */
    const GrainSizeDistribution* populationSizeDistribution(int f) const override;

    /** This function returns the dust mass \f$\mu_f\f$ per hydrogen atom for the fragment with
        index \f$f\f$ as it would be in the original dust mix being fragmented, without taking into
        account any additional fragment weights. */
    double populationMass(int f) const override;

    /** This function returns the total dust mass \f$\mu\f$ per hydrogen atom for all fragments
        combined in the original dust mix being fragmented, without taking into account any
        additional fragment weights. */
    double totalMass() const override;

    //=========== Exposing fragments (additional) ============

public:
    /** This function returns the equilibrium temperature of the dust population represented by the
        fragment with index \f$f\f$ when it would be embedded in a given radiation field. */
    double populationTemperature(int f, const Array& Jv) const;

    /** This function returns true if the human-readable identifier for the type of grain material
        represented by the fragment with index \f$f\f$ contains "Gra", "PAH" or "CM20", and false
        otherwise. */
    bool populationIsGraphite(int f) const;

    /** This function returns the average radius of a dust grain in the dust population represented
        by the fragment with index \f$f\f$. */
    double populationGrainRadius(int f) const;

    //======================== Data Members ========================

private:
    // list of fragments created in setupSelfAfter() and destructed automatically by being children
    vector<const DustMixFragment*> _fragments;
    int _numFrags{0};
};

////////////////////////////////////////////////////////////////////

#endif
