/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ABSORPTIONONLYMATERIALMIXDECORATOR_HPP
#define ABSORPTIONONLYMATERIALMIXDECORATOR_HPP

#include "MaterialMix.hpp"
#include "MultiGrainPopulationInterface.hpp"

////////////////////////////////////////////////////////////////////

/** AbsorptionOnlyMaterialMixDecorator is a material mix decorator that removes scattering from any
    material mix, passing through just the absorption properties. It can be used to evaluate the
    effects of scattering on observations. As a result, one can use SKIRT to test different
    versions of approximate radiation transfer:

    - a version that accounts for scattering out of the line of sight but neglects scattering into
    the line of sight; this can be obtained from the direct flux offered by the SKIRT instruments.

    - a version that neglects scattering both out of and into the line of sight; this can be
    obtained by removing all scattering from the material mix using this decorator.

    ## Fragmented dust mixes

    When working with fragmented dust mixes, first apply the FragmentDustMixDecorator to the
    original MultiGrainDustMix and then apply the AbsorptionOnlyMaterialMixDecorator to the
    FragmentDustMixDecorator. The other way around doesn't work because the
    FragmentDustMixDecorator explicitly expects a MultiGrainDustMix while the
    AbsorptionOnlyMaterialMixDecorator inherits directly from MaterialMix. This is not really a
    limitation as the intended goal is properly achieved.

    ## Limitations

    The current implementation of this decorator does not support secondary emission. To make this
    work, the decorator would need to pass through all material mix functions related to dust
    and/or gas emission, and the dynamic casting of material mix subclasses in the secondary
    emission framework would need to be replaced by the SimulationItem::interface() mechanism.

    Similarly, the current implementation of the decorator does not support dynamic dust
    destruction. To make this work, the decorator would need to implement the additional functions
    in the FragmentDustMixDecorator class provided for this purpose as a separate interface, and
    the DustDestructionRecipe class would need to use this interface rather than calling the
    FragmentDustMixDecorator class directly. */
class AbsorptionOnlyMaterialMixDecorator : public MaterialMix, public MultiGrainPopulationInterface
{
    ITEM_CONCRETE(AbsorptionOnlyMaterialMixDecorator, MaterialMix,
                  "a decorator that removes scattering from any material mix")
        ATTRIBUTE_TYPE_ALLOWED_IF(AbsorptionOnlyMaterialMixDecorator, "!Emission")
        ATTRIBUTE_TYPE_DISPLAYED_IF(AbsorptionOnlyMaterialMixDecorator, "Level3")

        PROPERTY_ITEM(materialMix, MaterialMix, "the material mix to be decorated")
        ATTRIBUTE_DEFAULT_VALUE(materialMix, "MeanInterstellarDustMix")

    ITEM_END()

    //======== Material type =======

public:
    /** This function returns the fundamental material type represented by the decorated material
        mix. */
    virtual MaterialType materialType() const override;

    //======== Capabilities =======

public:
    /** This function returns true if the decorated material mix supports polarization during
        scattering events, and false otherwise. */
    virtual bool hasPolarizedScattering() const override;

    /** This function returns true if the absorption of radiation for the decorated material mix is
        dichroic (i.e. the absorption cross section depends on the polarization state of incoming
        photon and the polarization state is adjusted during absorption), and false otherwise. */
    virtual bool hasPolarizedAbsorption() const override;

    /** This function returns true if scattering for the decorated material mix is resonant (such
        as for Lyman-alpha), and false otherwise. */
    virtual bool hasResonantScattering() const override;

    /** This function returns true if the extinction cross section (the sum of the absorption and
        scattering cross section) for the decorated material mix can be negative, and false
        otherwise. */
    virtual bool hasNegativeExtinction() const override;

    /** This function returns true if the cross sections returned by the decorated material mix may
        depend on the values of specific state variables other than the number density, and false
        otherwise. */
    virtual bool hasExtraSpecificState() const override;

    /** This function returns an enumeration indicating whether the decorated material mix offers
        an algorithm to update its specific medium state, and if so, when the updateSpecificState()
        function should be called. */
    virtual DynamicStateType hasDynamicMediumState() const override;

    //======== Medium state setup =======

public:
    /** This function returns the number and type of import parameters required by the decorated
        material mix as a list of SnapshotParameter objects. */
    virtual vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns a list of StateVariable objects describing the specific state
        variables used by the decorated material mix. */
    virtual vector<StateVariable> specificStateVariableInfo() const override;

    /** This function initializes any specific state variables requested by the decorated material
        mix through the specificStateVariableInfo() function except for the number density. */
    virtual void initializeSpecificState(MaterialState* state, double metallicity, double temperature,
                                         const Array& params) const override;

    //======== Medium state updates =======

public:
    /** If the decorated material mix has a dynamic medium state, i.e. if the
        hasDynamicMediumState() function returns anything other than \c None, this function is
        invoked for each spatial cell at the end of each relevant primary or secondary emission
        segment. Based on the specified radiation field, if needed, the function updates any values
        in the specific material state for this cell and medium component that may inform the local
        emission and/or extinction properties of the material. The function returns the update
        status as described for the UpdateStatus class. */
    virtual UpdateStatus updateSpecificState(MaterialState* state, const Array& Jv) const override;

    /** If the decorated material mix has a dynamic medium state, i.e. if the
        hasDynamicMediumState() function returns anything other than \c None, this function is
        invoked (once) \em after updateSpecificState() has been called for all spatial cells. The
        \em numCells, \em numUpdated and \em numNotConverged arguments specify respectively the
        number of spatial cells in the simulation, the number of cells updated during this update
        cycle, and the number of updated cells that have not yet converged. The \em
        currentAggregate and \em previousAggregate arguments provide the current and previous
        aggregate material states for this material mix (for more information, see the section on
        aggregation in the MediumState class header). Based on this information and any relevant
        user configuration options, the function returns true if the medium state is considered to
        be converged and false if not. */
    virtual bool isSpecificStateConverged(int numCells, int numUpdated, int numNotConverged,
                                          MaterialState* currentAggregate,
                                          MaterialState* previousAggregate) const override;

    //======== Low-level material properties =======

public:
    /** This function returns the mass per entity \f$\mu\f$ for the decorated material mix. */
    virtual double mass() const override;

    /** This function returns the default absorption cross section per entity
        \f$\varsigma^{\text{abs}}_{\lambda}\f$ at wavelength \f$\lambda\f$ for the decorated
        material mix. */
    virtual double sectionAbs(double lambda) const override;

    /** This function always returns zero as the default scattering cross section per entity
        \f$\varsigma^{\text{sca}}_{\lambda}\f$ at any wavelength \f$\lambda\f$, because this
        decorator blocks scattering. */
    virtual double sectionSca(double lambda) const override;

    /** This function returns the absorption cross section per entity of the decorated material mix
        as the default extinction cross section per entity, because this decorator blocks
        scattering: \f$\varsigma^{\text{ext}}_{\lambda} = \varsigma^{\text{abs}}_{\lambda}\f$ at
        all wavelengths \f$\lambda\f$. */
    virtual double sectionExt(double lambda) const override;

    //======== High-level photon life cycle =======

public:
    /** This function returns the absorption opacity \f$k^\text{abs}=n\varsigma^\text{abs}\f$ for
        the given wavelength, material state, and photon properties, for the decorated material
        mix. */
    virtual double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function always returns zero as the scattering opacity
        \f$k^\text{sca}=n\varsigma^\text{sca}\f$ for the given wavelength, material state, and
        photon properties, because this decorator blocks scattering. */
    virtual double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the absorption opacity for the decorated material mix as the
        extinction opacity, because this decorator blocks scattering:
        \f$k^\text{ext}=k^\text{abs}\f$ for the given wavelength, material state, and photon
        properties. */
    virtual double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    //============== Indicative temperature =============

public:
    /** This function returns an indicative temperature for the material represented by the
        specified material state and the decorated material mix, assuming an embedding radiation
        field specified by the mean intensities \f$(J_\lambda)_\ell\f$, if available. The
        interpretation of the indicative temperature depends heavily on the material type. */
    virtual double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //=========== Exposing multiple grain populations (MultiGrainPopulationInterface) ============

protected:
    /** This function is used by the interface() function to ensure that the receiving item can
        actually offer the specified interface. If the requested interface is the
        MultiGrainPopulationInterface, the implementation in this class returns true if the
        decorated material mix offers the MultiGrainPopulationInterface, and false otherwise. For
        other requested interfaces, the function invokes its counterpart in the base class. */
    bool offersInterface(const std::type_info& interfaceTypeInfo) const override;

public:
    /** If the decorated material mix offers the MultiGrainPopulationInterface, this function
        returns the result returned by the decorated material mix. Otherwise it throws a fatal
        error. */
    int numPopulations() const override;

    /** If the decorated material mix offers the MultiGrainPopulationInterface, this function
        returns the result returned by the decorated material mix. Otherwise it throws a fatal
        error. */
    string populationGrainType(int c) const override;

    /** If the decorated material mix offers the MultiGrainPopulationInterface, this function
        returns the result returned by the decorated material mix. Otherwise it throws a fatal
        error. */
    double populationBulkDensity(int c) const override;

    /** If the decorated material mix offers the MultiGrainPopulationInterface, this function
        returns the result returned by the decorated material mix. Otherwise it throws a fatal
        error. */
    Range populationSizeRange(int c) const override;

    /** If the decorated material mix offers the MultiGrainPopulationInterface, this function
        returns the result returned by the decorated material mix. Otherwise it throws a fatal
        error. */
    const GrainSizeDistribution* populationSizeDistribution(int c) const override;

    /** If the decorated material mix offers the MultiGrainPopulationInterface, this function
        returns the result returned by the decorated material mix. Otherwise it throws a fatal
        error. */
    double populationMass(int c) const override;

    /** If the decorated material mix offers the MultiGrainPopulationInterface, this function
        returns the result returned by the decorated material mix. Otherwise it throws a fatal
        error. */
    double totalMass() const override;
};

////////////////////////////////////////////////////////////////////

#endif
