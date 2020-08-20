/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FRAGMENTDUSTMIXDECORATOR_HPP
#define FRAGMENTDUSTMIXDECORATOR_HPP

#include "MultiGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** TO DO. */
class FragmentDustMixDecorator : public MaterialMix, public MultiGrainPopulationInterface
{
    ITEM_CONCRETE(FragmentDustMixDecorator, MaterialMix,
                  "a dust mix decorator that manages separate densities for fragments of another dust mix")
        ATTRIBUTE_TYPE_INSERT(FragmentDustMixDecorator, "Dust")
        ATTRIBUTE_TYPE_DISPLAYED_IF(FragmentDustMixDecorator, "Level2")

        PROPERTY_ITEM(dustMix, MultiGrainDustMix, "the dust mix to be fragmented")
        ATTRIBUTE_DEFAULT_VALUE(dustMix, "ThemisDustMix")

        PROPERTY_BOOL(fragmentSizeBins, "fragment each dust population into its grain size bins")
        ATTRIBUTE_DEFAULT_VALUE(fragmentSizeBins, "false")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** TO DO. */
    void setupSelfAfter() override;

    //======== Material type =======

public:
    /** This function returns the fundamental material type represented by this material mix. For
        this dust mix decorator, it returns MaterialType::Dust. */
    MaterialType materialType() const override;

    //======================== Capabilities =======================

public:
    /** This function returns a flag indicating whether the material mix supports polarization
        during scattering events or not. For this dust mix decorator, the function returns the
        value obtained from the corresponding function for the dust mix being fragmented. */
    bool hasPolarizedScattering() const override;

    /** This function returns a flag indicating whether this material mix supports stochastic
        heating of dust grains. For this dust mix decorator, the function returns true because
        all dust mixes being fragmented support stochastic heating by definition. */
    bool hasStochasticDustEmission() const override;

    //======== Medium state setup =======

public:
    /** This function returns the number and type of import parameters required by this particular
        material mix as a list of SnapshotParameter objects. TO DO. */
    vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns a list of StateVariable objects describing the specific state
        variables used by the receiving material mix. TO DO. */
    vector<StateVariable> specificStateVariableInfo() const override;

    /** This function initializes any specific state variables requested by this material mix
        through the specificStateVariableInfo() function except for the number density. TO DO. */
    void initializeSpecificState(MaterialState* state, double temperature, const Array& params) const override;

    //======== Low-level material properties =======

public:
    /** This function returns the dust mass \f$\mu\f$ per hydrogen atom for this dust mix. TO DO.
        */
    double mass() const override;

    /** This function returns the absorption cross section per entity
        \f$\varsigma^{\text{abs}}_{\lambda}\f$ at wavelength \f$\lambda\f$. TO DO. */
    double sectionAbs(double lambda) const override;

    /** This function returns the scattering cross section per entity
        \f$\varsigma^{\text{sca}}_{\lambda}\f$ at wavelength \f$\lambda\f$. TO DO. */
    double sectionSca(double lambda) const override;

    /** This function returns the total extinction cross section per entity
        \f$\varsigma^{\text{ext}}_{\lambda} = \varsigma^{\text{abs}}_{\lambda} +
        \varsigma^{\text{sca}}_{\lambda}\f$ at wavelength \f$\lambda\f$. TO DO. */
    double sectionExt(double lambda) const override;

    /** This function returns the scattering asymmetry parameter \f$g_\lambda =
        \left<\cos\theta\right>\f$ at wavelength \f$\lambda\f$, which is used with the
        HenyeyGreenstein scattering mode. TO DO. */
    double asymmpar(double lambda) const override;

    //======== High-level photon life cycle =======

    /** This function returns the absorption opacity \f$k^\text{abs}=n\varsigma^\text{abs}\f$ for
        the given wavelength and material state. The photon properties are not used. TO DO. */
    double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the scattering opacity \f$k^\text{sca}=n\varsigma^\text{sca}\f$ for
        the given wavelength and material state. The photon properties are not used. TO DO. */
    double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the extinction opacity \f$k^\text{ext}=k^\text{abs}+k^\text{sca}\f$
        for the given wavelength and material state. The photon properties are not used. TO DO. */
    double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function calculates the contribution of the medium component associated with this
        material mix to the peel-off photon luminosity, polarization state, and wavelength shift
        for the given wavelength, geometry, material state, and photon properties. TO DO. */
    void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, double w, Direction bfkobs,
                           Direction bfky, const MaterialState* state, PhotonPacket* pp) const override;

    /** This function performs a scattering event on the specified photon packet in the spatial
        cell and medium component represented by the specified material state and the receiving
        material mix. TO DO. */
    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    /** This function returns an indicative temperature of the material mix when it would be
        embedded in a given radiation field. TO DO. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //======== Secondary emission =======

public:
    /** This function returns the emissivity spectrum per hydrogen atom \f$\varepsilon_{\ell'}\f$
        of the dust mix (or rather of the representative grain population corresponding to the dust
        mix) when it would be embedded in a given radiation field, assuming that the dust grains
        are in local thermal equilibrium. The input and output arrays are discretized on the
        wavelength grids returned by the Configuration::radiationFieldWLG() and
        Configuration::dustEmissionWLG() functions, repectively. TO DO. */
    Array emissivity(const Array& Jv) const override;

    //=========== Exposing multiple grain populations (MultiGrainPopulationInterface) ============

public:
    /** This function returns the number of dust grain populations (with indices \f$c\f$) in this
        dust mix. TO DO. */
    int numPopulations() const override;

    /** This function returns a brief human-readable identifier for the type of grain material
        represented by the population with index \f$c\f$. TO DO. */
    string populationGrainType(int c) const override;

    /** This function returns the minimum and maximum grain sizes \f$a_{\text{min},c},
        a_{\text{max},c}\f$ for the population with index \f$c\f$. TO DO. */
    Range populationSizeRange(int c) const override;

    /** This function returns the grain size distribution object for the population with index
        \f$c\f$. TO DO. */
    const GrainSizeDistribution* populationSizeDistribution(int c) const override;

    /** This function returns the dust mass \f$\mu_c\f$ per hydrogen atom for the population with
        index \f$c\f$. TO DO. */
    double populationMass(int c) const override;

    /** This function returns the total dust mass \f$\mu_c\f$ per hydrogen atom for all populations
        combined. TO DO. */
    double totalMass() const override;

    //======================== Data Members ========================

private:
    // all data members are precalculated in setupSelfAfter()
};

////////////////////////////////////////////////////////////////////

#endif
