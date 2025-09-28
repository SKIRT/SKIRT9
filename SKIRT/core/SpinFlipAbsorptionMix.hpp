/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPINFLIPABSORPTIONMIX_HPP
#define SPINFLIPABSORPTIONMIX_HPP

#include "MaterialMix.hpp"

////////////////////////////////////////////////////////////////////

/** The SpinFlipAbsorptionMix class describes the material properties related to the 21 cm
    spin-flip transition in neutral atomic hydrogen, including only absorption (i.e. omitting
    emission processes). The 21 cm opacity in a given cell is determined from gas properties
    defined in the input model (gas temperature).

    <b>Configuring the simulation</b>

    Simulations of 21 cm the spin-flip transition usually include 21 cm emission (see the
    SpinFlipSEDFamily class) in addition to a medium component configured with the spin flip
    material mix (this class). Hence, the simulation mode should be set to "ExtinctionOnly".

    As the opacity of the 21 cm line depends on the gas temperature as well as the neutral atomic
    hydrogen density \f$n_\mathrm{HI}\f$, this information is read from an input file by
    associating the SpinFlipAbsorptionMix with a subclass of ImportedMedium. For that medium
    component, the ski file attribute \em importTemperature <b>must</b> be set to 'true', and \em
    importMetallicity and \em importVariableMixParams must be left at 'false'. For example, if bulk
    velocities are also imported for this medium component (i.e. \em importVelocity is 'true'), the
    column order would be \f[ ..., T, v_\mathrm{x}, v_\mathrm{y}, v_\mathrm{z} \f]

    For basic testing purposes, the SpinFlipAbsorptionMix can also be associated with a geometric
    medium component. The geometry then defines the spatial density distribution (i.e.
    \f$n_\mathrm{HI}\f$), and this class offers a configuration property to specify a fixed default
    temperature value that will be used across the spatial domain.

    <b>Absorption</b>

    Following Draine 2011 Chapter 8, and substituting wavelengths for frequencies, the
    monochromatic absorption cross section per neutral hydrogen atom \f$\varsigma(\lambda)\f$ at
    frequency \f$\lambda\f$ can be written as \f[\varsigma^\text{abs}(\lambda) =
    \frac{3}{32\pi}\,A_\mathrm{SF}\,\frac{hc \,\lambda_\mathrm{SF}} {k_\mathrm{B}T_\mathrm{s}}
    \times \frac{\lambda_\mathrm{SF}} {\sqrt{2\pi}\,\sigma} \,\exp\left(-\frac{u^2(\lambda)}
    {2\sigma^2}\right) \f] where \f$T_\mathrm{s}\f$ is the spin temperature of the hydrogen gas
    (see below), \f$\sigma=\sqrt{k_B T/m_\mathrm{p}}\f$ is the velocity dispersion corresponding to
    the thermal motion of the atoms, and \f$u(\lambda)=c(\lambda-\lambda_\mathrm{SF})/\lambda\f$ is
    the frequency deviation from the line center in velocity units.

    According to Kim, Ostriker and Kim 2014 (ApJ 786:64), the spin temperature \f$T_\mathrm{s}\f$
    is essentially equal to the kinetic gas temperature \f$T\f$ for low gas temperatures
    \f$T<1000~\mathrm{K}\f$. For higher gas temperatures, it also depends on other gas properties,
    with an upper bound that never exceeds a fixed maximum value, i.e. \f$T_\mathrm{s,upper}
    \lesssim 5000~\mathrm{K}\f$ (see their Figure 2).

    To avoid the need for defining additional gas properties in the input model, we here use an
    approximate upper bound for the spin temperature given by \f[ T_\mathrm{s}= 6000~\mathrm{K}
    \times \left( 1-\mathrm{e}^{-T/5000~\mathrm{K} } \right).\f] Because the absorption cross
    section is inversely proportional to the spin temperature, this yields an approximate lower
    bound for the cross section.

    */
class SpinFlipAbsorptionMix : public MaterialMix
{
    ITEM_CONCRETE(SpinFlipAbsorptionMix, MaterialMix, "A gas mix supporting the spin-flip 21 cm hydrogen absorption")

        PROPERTY_DOUBLE(defaultTemperature, "the default temperature of the gas")
        ATTRIBUTE_QUANTITY(defaultTemperature, "temperature")
        ATTRIBUTE_MIN_VALUE(defaultTemperature, "[3")  // gas temperature must be above local Universe T_CMB
        ATTRIBUTE_MAX_VALUE(defaultTemperature, "1e9]")
        ATTRIBUTE_DEFAULT_VALUE(defaultTemperature, "1e4")
        ATTRIBUTE_DISPLAYED_IF(defaultTemperature, "Level2")

    ITEM_END()

    //======== Capabilities =======

public:
    /** This function returns the fundamental material type represented by this material mix, which
        is MaterialType::Gas. */
    MaterialType materialType() const override;

    /** This function returns true, indicating that the cross sections returned by this material
        mix depend on the values of specific state variables other than the number density. */
    bool hasExtraSpecificState() const override;

    //======== Medium state setup =======

public:
    /** This function returns a list of StateVariable objects describing the specific state
        variables used by the receiving material mix. For this class, the function returns a list
        containing descriptors for the properties defined in the input model, namely number density
        and temperature. */
    vector<StateVariable> specificStateVariableInfo() const override;

    /** This function initializes the specific state variables requested by this material mix
        through the specificStateVariableInfo() function except for the number density. For this
        class, the function initializes the temperature to the specified imported value, or if not
        available, to the user-configured default value. */
    void initializeSpecificState(MaterialState* state, double metallicity, double temperature,
                                 const Array& params) const override;

    //======== Low-level material properties =======

public:
    /** This function returns the mass of a neutral hydrogen atom. */
    double mass() const override;

    /** This function returns the 21 cm absorption cross section per neutral hydrogen atom at the
        given wavelength and using the default gas properties configured for this material mix. */
    double sectionAbs(double lambda) const override;

    /** This function returns the 21 scattering cross section per neutral hydrogen atom, which is
        trivially zero for all wavelengths. */
    double sectionSca(double lambda) const override;

    /** This function returns the total 21 cm extinction cross section per neutral hydrogen atom at
        the given wavelength and using the default gas properties configured for this material mix.
        The extinction cross section is identical to the absorption cross section because the
        scattering cross section is zero. */
    double sectionExt(double lambda) const override;

    //======== High-level photon life cycle =======

    /** This function returns the 21 cm absorption opacity \f$k^\text{abs}= n_mathrm{HI}
        \varsigma^\text{abs}\f$ for the given wavelength and material state. The photon packet
        properties are not used. */
    double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the 21 cm scattering opacity \f$k^\text{sca}\f$ which is trivially
        zero at all wavelengths. */
    double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the 21 cm extinction opacity \f$k^\text{ext}=k^\text{abs}\f$ for the
        given wavelength and material state, which equals the absorption opacity because the
        scattering opacity is zero. The photon packet properties are not used. */
    double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    //======== Temperature =======

    /** This function returns an indicative temperature of the material mix when it would be
        embedded in a given radiation field. The implementation in this class ignores the radiation
        field and returns the temperature stored in the specific state for the relevant spatial
        cell and medium component. Because the hydrogen temperature is not calculated
        self-consistently in our treatment, this value corresponds to the temperature defined by
        the input model at the start of the simulation. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;
};

////////////////////////////////////////////////////////////////////

#endif
