/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPINFLIPHYDROGENGASMIX_HPP
#define SPINFLIPHYDROGENGASMIX_HPP

#include "EmittingGasMix.hpp"

////////////////////////////////////////////////////////////////////

/** The SpinFlipHydrogenGasMix class describes the material properties related to the 21 cm
    spin-flip transition in neutral atomic hydrogen, including emission and absorption. The 21 cm
    emission luminosity and self-absorption opacity in a given cell are determined from gas
    properties defined in the input model (total hydrogen number density, gas metallity, gas
    temperature, neutral hydrogen fraction) and the local UV radiation field calculated by the
    simulation (taking into account dust extinction).

    <b>Partitioning</b>

    While the input model is expected to define the neutral (i.e. non-ionized) hydrogen fraction,
    the partitioning of the hydrogen gas into its atomic and molecular components must be
    determined based on the the UV radiation field calculated by the simulation. Since UV radiation
    in the Lyman-Werner band (912-1110 Angstrom) dissociates molecular hydrogen, hydrogen
    partitioning schemes usually have the UV field in that wavelength range as an important
    parameter. For this purpose, this class samples the radiation field
    \f$J_\lambda(\lambda_\mathrm{UV})\f$ at \f$\lambda_\mathrm{UV}=1000\,\mathrm{Angstrom}\f$ The
    wavelength bin index for \f$\lambda_\mathrm{UV}\f$ in the radiation field grid is determined
    during setup.

    Currently, this class implements the partioning scheme described by Gnedin and Kravtsov 2011
    (ApJ 728:88). This scheme estimates the mass fraction of molecular hydrogen compared to the
    total hydrogen mass, i.e. \f$f_\mathrm{H2} = M_\mathrm{H2}/M_H\f$, where \f$M_\mathrm{H} =
    M_\mathrm{HI} + M_\mathrm{H2} + M_\mathrm{HII}\f$ denotes the total (atomic, molecular, and
    ionized) hydrogen mass. The scheme has the following three input parameters:

    - The total hydrogen number density \f$n_\mathrm{H}\f$ defined in the input model.

    - The dust-to-gas ratio \f$D\f$ relative to the representative Galactic value. Because the
    dust-to-gas ratio can be assumed to scale with metallicity, this is roughly equivalent to the
    gas metallicity relative to the solar value. In other words, \f$D=Z/0.0127\f$, where \f$Z\f$ is
    defined in the input model.

    - The dimensionless UV radiation field \f$U = J_\lambda(\lambda_\mathrm{UV}) /
    J^\mathrm{MW}_\lambda(\lambda_\mathrm{UV})\f$, obtained by normalizing the sampled UV radiation
    field to the representative Galactic value, which is taken to be
    \f$J^\mathrm{MW}_{E}(\lambda_\mathrm{UV}) = 10^{6} \,\mathrm{photons} \,\mathrm{cm}^{-2}
    \,\mathrm{s}^{-1} \,\mathrm{sr}^{-1} \,\mathrm{eV}^{-1}\f$. The latter value is converted at
    compile time from the specified photons-per-energy units to the internal SKIRT
    energy-per-wavelength units.

    The Gnedin and Kravtsov 2011 partitioning scheme is then defined by their equations (6) and
    following, which are replicated below using our notation.

    \f[ f_\mathrm{H2} = \frac{1}{1 + \exp(-4x-3x^3)} \f]

    \f[ x = \Lambda^{3/7}\,\ln\left( \frac{D\,n_\mathrm{H}}{\Lambda\,n_*}\right) \f]

    \f[ \Lambda = \ln\left( 1 + g D^{3/7} (U/15)^{4/7} \right) \f]

    \f[ g = \frac{1+\alpha s + s^2}{1+s} \f]

    \f[ s = \frac{0.04}{D_* + D} \f]

    \f[ \alpha = \frac{5\,U/2}{1+(U/2)^2} \f]

    \f[ n_* = 25~\mathrm{cm}^{-3} \f]

    \f[ D_* = 1.5 \times 10^{-3} \times \ln\left(1 + (3 U)^{1.7} \right) \f]

    Given the total hydrogen number density \f$n_\mathrm{H}\f$ and the neutral hydrogen fraction
    \f$f_\mathrm{HI+H2} = (M_\mathrm{HI} + M_\mathrm{H2})/M_H\f$ defined in the input model, the
    atomic number density can then be easily derived as \f$n_\mathrm{HI} = n_\mathrm{H}
    (f_\mathrm{HI+H2} - f_\mathrm{H2})\f$.

    <b>Emission</b>

    Following Draine 2011 Chapter 8, the integrated luminosity \f$L\f$ of the 21 cm line for a
    given spatial cell can be written as \f[ L = \frac{3}{4} \,A_\mathrm{SF}
    \,\frac{hc}{\lambda_\mathrm{SF}} \,n_\mathrm{HI}V \f] where \f$A_\mathrm{SF}=2.8843 \times
    10^{-15} \,s^{-1}\f$ is the Einstein coefficient of the 21cm spin-flip transition,
    \f$\lambda_\mathrm{SF} = 21.10611405413\times 10^{-2}\,\mathrm{m}\f$ is its wavelength,
    \f$n_\mathrm{HI}\f$ is the number density of atomic hydrogen in the cell, and \f$V\f$ is the
    cell volume.

    The 21cm line is extremely narrow so that the emerging line profile is dominated by the Doppler
    shift caused by the thermal motion of the atoms. The gas temperature defined in the input model
    is used to randomly shift the wavelength of each emitted photon packet according to a
    Maxwell-Boltzmann distribution. This happens outside of this class in the secondary source
    machinery.

    <b>Absorption</b>

    Again following Draine 2011 Chapter 8, the monochromatic absorption cross section
    \f$\varsigma(\lambda)\f$ at wavelength \f$\lambda\f$ can be written as \f[\varsigma(\lambda) =
    \frac{3}{32\pi}\,\frac{A_\mathrm{SF}\,hc \,\lambda_\mathrm{SF}} {k_\mathrm{B}T_\mathrm{s}}
    \,\phi(\lambda)\f] where \f$T_\mathrm{s}\f$ is the spin temperature of the hydrogen gas (see
    below) and \f$\phi(\lambda)\f$ is a Gaussian distribution representing the Doppler shift caused
    by the thermal motion of the atoms, with \f$\int\phi(\lambda) \,\mathrm{d}\lambda=1\f$.

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

    <b>Configuring the simulation</b>

    simulation mode

    RF bins

    default values

    The spatial distributions for the gas density, metallicity, and temperature and for the local
    dust-to-gas ratio must be defined by the input model and are considered to be constant during
    the simulation. In this context, this material mix offers configuration properties to specify
    default values for these quantities that will be used by geometric media across the spatial
    domain. */
class SpinFlipHydrogenGasMix : public EmittingGasMix
{
    ITEM_CONCRETE(SpinFlipHydrogenGasMix, EmittingGasMix,
                  "A gas mix supporting the spin-flip 21 cm hydrogen transition")
        ATTRIBUTE_TYPE_INSERT(SpinFlipHydrogenGasMix, "CustomMediumState")

        PROPERTY_DOUBLE(defaultMetallicity, "the default metallicity of the gas")
        ATTRIBUTE_MIN_VALUE(defaultMetallicity, "]0")
        ATTRIBUTE_MAX_VALUE(defaultMetallicity, "0.2]")
        ATTRIBUTE_DEFAULT_VALUE(defaultMetallicity, "0.02")
        ATTRIBUTE_DISPLAYED_IF(defaultMetallicity, "Level2")

        PROPERTY_DOUBLE(defaultTemperature, "the default temperature of the gas")
        ATTRIBUTE_QUANTITY(defaultTemperature, "temperature")
        ATTRIBUTE_MIN_VALUE(defaultTemperature, "[3")  // gas temperature must be above local Universe T_CMB
        ATTRIBUTE_MAX_VALUE(defaultTemperature, "1e9]")
        ATTRIBUTE_DEFAULT_VALUE(defaultTemperature, "1e4")
        ATTRIBUTE_DISPLAYED_IF(defaultTemperature, "Level2")

        PROPERTY_DOUBLE(defaultDustToGasRatio, "the default dust-to-gas ratio")
        ATTRIBUTE_MIN_VALUE(defaultDustToGasRatio, "]0")
        ATTRIBUTE_MAX_VALUE(defaultDustToGasRatio, "0.2]")
        ATTRIBUTE_DEFAULT_VALUE(defaultDustToGasRatio, "0.01")
        ATTRIBUTE_DISPLAYED_IF(defaultDustToGasRatio, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function determines the radiation field wavelength bin containing the UV field
        strength. */
    void setupSelfBefore() override;

    //======== Capabilities =======

public:
    /** This function returns the fundamental material type represented by this material mix, which
        is MaterialType::Gas. */
    MaterialType materialType() const override;

    /** This function returns true, indicating that the cross sections returned by this material
        mix depend on the values of specific state variables other than the number density. */
    bool hasExtraSpecificState() const override;

    /** This function returns true, indicating that this material has a semi-dynamic medium state.
        */
    bool hasSemiDynamicMediumState() const override;

    /** This function returns true, indicating that this material supports secondary line emission
        from gas. */
    bool hasLineEmission() const override;

    //======== Medium state setup =======

public:
    /** This function returns the number and type of import parameters required by this particular
        material mix as a list of SnapshotParameter objects. For this class, the function returns a
        descriptor for the dust-to-gas ratio import parameter. Importing metallicity and
        temperature should be enabled through the corresponding standard configuration flags. */
    vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns a list of StateVariable objects describing the specific state
        variables used by the receiving material mix. For this class, the function returns a list
        containing descriptors for number density, metallicity, temperature, a custom variable to
        hold the dust-to-gas-ratio, and a custom variable to hold the UV field strength derived
        from the radiation field when the semi-dynamic medium state is updated. */
    vector<StateVariable> specificStateVariableInfo() const override;

    /** This function initializes the specific state variables requested by this fragmented dust
        mix through the specificStateVariableInfo() function except for the number density. For
        this class, the function initializes the temperature, metallicity and dust-to-gas ratio to
        the specified imported values, or if not available, to the user-configured default values.
        The UV field strength is set to zero. */
    void initializeSpecificState(MaterialState* state, double metallicity, double temperature,
                                 const Array& params) const override;

    //======== Medium state updates =======

    /** Based on the specified radiation field, the function obtains the UV field strength and
        stores it in the medium state for the specified cell. The function returns true if the
        medium state has indeed be changed, and false otherwise. */
    bool updateSpecificState(MaterialState* state, const Array& Jv) const override;

    //======== Low-level material properties =======

public:
    /** This function returns the mass of a neutral hydrogen atom. */
    double mass() const override;

    /** This function returns the 21 cm absorption cross section per hydrogen atom at the given
        wavelength and using the default gas properties configured for this material mix. */
    double sectionAbs(double lambda) const override;

    /** This function returns the 21 scattering cross section per hydrogen atom, which is trivially
        zero for all wavelengths. */
    double sectionSca(double lambda) const override;

    /** This function returns the total 21 cm extinction cross section per hydrogen atom at the
        given wavelength and using the default gas properties configured for this material mix. The
        extinction cross section is identical to the absorption cross section because the
        scattering cross section is zero. */
    double sectionExt(double lambda) const override;

    //======== High-level photon life cycle =======

    /** This function returns the 21 cm absorption opacity \f$k^\text{abs}=n\varsigma^\text{abs}\f$
        for the given wavelength and material state. The photon packet properties are not used. */
    double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the 21 cm scattering opacity \f$k^\text{sca}=n\varsigma^\text{sca}\f$
        which is trivially zero at all wavelengths. */
    double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the 21 cm extinction opacity
        \f$k^\text{ext}=k^\text{abs}+k^\text{sca}\f$ for the given wavelength and material state.
        The photon packet properties are not used. The extinction opacity equals the absorption
        opacity. */
    double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function throws a fatal error because the 21 cm line does not scatter. */
    void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, double w, Direction bfkobs,
                           Direction bfky, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function throws a fatal error because the 21 cm line does not scatter. */
    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======== Secondary emission =======

    /** This function returns a list including a single item: the line center of the 21 cm hydrogen
        spinflip transition. */
    Array lineEmissionCenters() const override;

    /** This function returns a list including a single item: the mass of the particle emitting the
        21 cm line, i.e. the hydrogen atom. */
    Array lineEmissionMasses() const override;

    /** This function returns  a list including a single item: the 21 cm line luminosity
        in the spatial cell and medium component represented by the specified material state and
        the receiving material mix when it would be embedded in the specified radiation field. */
    Array lineEmissionSpectrum(const MaterialState* state, const Array& Jv) const override;

    //======== Temperature =======

    /** This function returns an indicative temperature of the material mix when it would be
        embedded in a given radiation field. The implementation in this class ignores the radiation
        field and returns the temperature stored in the specific state for the relevant spatial
        cell and medium component. Because the hydrogen temperature is not calculated
        self-consistently in our treatment, this value corresponds to the temperature defined by
        the input model at the start of the simulation. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //======================== Data Members ========================

private:
    int _indexUV{-1};  // index in simulation's RF WLG for the bin containing 1000 A
};

////////////////////////////////////////////////////////////////////

#endif
