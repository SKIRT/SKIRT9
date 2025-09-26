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
    properties defined in the input model (neutral hydrogen volume density, gas metallity, gas
    temperature, and neutral hydrogen surface density) and the local UV radiation field calculated
    by the simulation taking into account dust extinction.

    <b>Configuring the simulation</b>

    Simulations of the 21 cm the spin-flip transition usually include primary sources and a dust
    medium in addition to a medium component configured with the spin flip material mix (this
    class). During primary emission, the dust medium determines the UV radiation field, allowing
    the spin flip material mix to calculate the 21 cm line luminosity and absorption cross section
    for use during secondary emission. This calculation happens in the updateSpecificState()
    function, which is invoked at the end of primary emission (because the material mix advertises
    that it has a secondary dynamic medium state). This does imply that the 21 cm line absorption
    remains zero during primary emission, which is not a problem as long as the primary sources
    don't emit in the radio wavelength range. Thus, unless dust opacities are very high, there is
    no need for iteration over primary nor secondary emission.

    If one is not interested in dust emission (i.e. only in the 21 cm line emission), the
    simulation mode can be set to "GasEmission" and the radiation field wavelength grid can be
    limited to a single bin in the UV wavelength range. This bin should be configured to cover the
    full Lyman-Werner band (912-1110 Å) so that it accurately captures all radiation that
    dissociates molecular hydrogen. Similarly, instrument wavelength grids can be limited to a
    relevant range around the 21 cm line center.

    To calculate both the dust emission spectrum and the 21 cm line emission in the same
    simulation, the simulation mode must be set to "DustAndGasEmission" and the radiation field
    wavelength grid must now include and properly resolve the UV, optical, and infrared wavelength
    range. Separate instruments can be configured for the relevant wavelength ranges, e.g. using a
    logarithmic grid for the continuum spectrum and a linear grid for the 21 cm line profile.

    The input model must provide values for the spatial distribution of four medium properties that
    remain constant during the simulation:

    - The neutral hydrogen volume number density \f$n_\mathrm{HI+H2} = n_\mathrm{HI} +
    2\,n_\mathrm{H2}\f$, i.e. including non-ionized atomic and molecular hydrogen. Note the factor
    2 in this equation; \f$n_\mathrm{HI+H2}\f$ specifies the number of protons rather than the
    number of hydrogen-like particles.

    - The metallicity \f$Z\f$ of the gas.

    - The kinetic temperature \f$T\f$ of the gas.

    - The neutral hydrogen surface mass density \f$\Sigma_\mathrm{HI+H2}\f$, required by the scheme
    that partitions the neutral hydrogen into atomic and molecular fractions (see below).

    Most often, this information will be read from an input file by associating the spin flip
    material mix with a subclass of ImportedMedium. For that medium component, the ski file
    attributes \em importMetallicity and \em importTemperature <b>must</b> be set to 'true', and
    \em importVariableMixParams must be left at 'false'. The additional column required by the spin
    flip material mix (neutral hydrogen surface mass density) is automatically imported and is
    expected <b>after</b> all other columns. For example, if bulk velocities are also imported for
    this medium component (i.e. \em importVelocity is 'true'), the column order would be \f[ ...,
    n_\mathrm{HI+H2}, Z, T, v_\mathrm{x}, v_\mathrm{y}, v_\mathrm{z}, \Sigma_\mathrm{HI+H2} \f]

    Depending on the specific ImportedMedium subclass being used, the neutral hydrogen density
    distribution in the first column can be specified as a mass or number density, or as a
    volume-integrated mass or number for each cell or particle. In all cases, the value is
    automatically converted to the appropriate volume number density expected by this class. In
    contrast, the neutral hydrogen surface density must always be specified as a surface \em mass
    density, with default units of Msun/pc2.

    \note It is customary to obtain \f$\Sigma_\mathrm{HI+H2}\f$ by multiplying
    \f$n_\mathrm{HI+H2}\f$ with a specific length scale (and converting from number to mass
    density). Some authors use the Sobolev length, but most use the Jeans length. The Jeans length
    can be easily computed from the gas cell data (gas mass density, internal energy, hydrogen mass
    fraction, and electron fraction). This rough approximation broadly matches the interpretation
    of surface density in the partitioning scheme used by this class (see below).

    For basic testing purposes, the spin flip material mix can also be associated with a geometric
    medium component. The geometry then defines the spatial neutral hydrogen density distribution
    (i.e. \f$n_\mathrm{HI+H2}\f$), and the spin flip material mix offers configuration properties
    to specify fixed default values for the other properties that will be used across the spatial
    domain.

    <b>Partitioning scheme</b>

    While the input model defines the neutral (i.e. non-ionized) hydrogen density distribution, the
    partitioning of the hydrogen gas into its atomic and molecular components must be determined
    based on the the UV radiation field calculated by the simulation. Since UV radiation in the
    Lyman-Werner band (912-1110 Å) dissociates molecular hydrogen, a hydrogen partitioning scheme
    usually has the UV field in that wavelength range as an important parameter. For this purpose,
    this class samples the radiation field \f$J_\lambda(\lambda_\mathrm{UV})\f$ at
    \f$\lambda_\mathrm{UV}=1000\,\mathrm{Å}\f$. The wavelength bin index for
    \f$\lambda_\mathrm{UV}\f$ in the radiation field grid is determined during setup.

    This class implements the column-density based partioning scheme of Gnedin & Draine 2014 (Apj,
    795, 37), which is an update of Gnedin & Kravtsov 2011 (ApJ, 728, 88). The mechanism for
    estimating the length scale of a gas cell is taken from Diemer et al. 2018 (ApJS, 238, 33).

    The scheme has the following input parameters:

    - The neutral hydrogen surface mass density \f$\Sigma_\mathrm{HI+H2}\f$ defined in the input
    model.

    - The dust-to-gas ratio \f$D_\mathrm{MW}\f$ relative to the representative Milky Way value.
    Because the dust-to-gas ratio can be assumed to scale with metallicity, this is roughly
    equivalent to the gas metallicity relative to the solar value. In other words,
    \f$D_\mathrm{MW}=Z/0.0127\f$, where \f$Z\f$ is the metallicity defined in the input model.

    - The dimensionless UV radiation field \f$U_\mathrm{MW} = J_\lambda(\lambda_\mathrm{UV}) /
    J^\mathrm{MW}_\lambda(\lambda_\mathrm{UV})\f$, obtained by normalizing the sampled UV radiation
    field to the representative Milky Way value, which is taken to be
    \f$J^\mathrm{MW}_{E}(\lambda_\mathrm{UV}) = 10^{6} \,\mathrm{photons} \,\mathrm{cm}^{-2}
    \,\mathrm{s}^{-1} \,\mathrm{sr}^{-1} \,\mathrm{eV}^{-1}\f$. The latter value is converted at
    compile time from the specified photons-per-energy cgs units to the internal SKIRT
    energy-per-wavelength SI units using \f[ J^\mathrm{MW}_\lambda(\lambda_\mathrm{UV}) = 10^4 \,
    \frac{(hc)^2}{q_\mathrm{el}\lambda^3_\mathrm{UV}} J^\mathrm{MW}_{E}(\lambda_\mathrm{UV}). \f]

    - The length scale \f$L_\mathrm{cell}\f$, introduced to obtain converged results when applying
    the partitioning scheme to simulations of varying resolutions. An estimate for this value is
    obtained using \f$L_\mathrm{cell}=V_\mathrm{cell}^{1/3}\f$, where \f$V_\mathrm{cell}\f$ denotes
    the volume of the spatial grid cell under consideration.

    With these input parameters, the partitioning scheme can be summarized as follows:

    \f[ S=L_\mathrm{cell}/100\,\mathrm{pc} \f]

    \f[ D_\ast=0.17\frac{2+S^5}{1+S^5} \f]

    \f[ g=\sqrt{D_\mathrm{MW}^2+D_\ast^2} \f]

    \f[ \Sigma_c=5\times10^7\,\mathrm{M}_\odot\mathrm{kpc}^{-2}\times
        \Bigl(\frac{\sqrt{0.001+0.1U_\mathrm{MW}}}{g(1+1.69\sqrt{0.001+0.1U_\mathrm{MW}}}\Bigr) \f]

    \f[ \alpha=0.5+\frac{1}{1+\sqrt{U_\mathrm{MW}D_\mathrm{MW}^2/600}} \f]

    \f[ R_\mathrm{mol} = \Bigl(\frac{\Sigma_\mathrm{HI+H_2}}{\Sigma_c}\Bigr)^\alpha \f]

    \f[ f_\mathrm{mol} = \frac{R_\mathrm{mol}}{R_\mathrm{mol}+1} \f]

    where \f$f_\mathrm{mol}\f$ is the fraction of molecular hydrogen in the neutral hydrogen gas,
    and \f$f_\mathrm{atom} = 1 - f_\mathrm{mol}\f$ is the corresponding atomic hydrogen fraction.

    <b>Emission</b>

    Following Draine 2011 Chapter 8, the integrated luminosity \f$L\f$ of the 21 cm line for a
    given spatial cell can be written as \f[ L = \frac{3}{4} \,A_\mathrm{SF}
    \,\frac{hc}{\lambda_\mathrm{SF}} \,n_\mathrm{HI}V \f] where \f$A_\mathrm{SF}=2.8843 \times
    10^{-15} \,s^{-1}\f$ is the Einstein coefficient of the 21cm spin-flip transition,
    \f$\lambda_\mathrm{SF} = 21.10611405413\times 10^{-2}\,\mathrm{m}\f$ is its wavelength,
    \f$n_\mathrm{HI}\f$ is the number density of atomic hydrogen in the cell, and \f$V\f$ is the
    cell volume.

    The 21 cm line is extremely narrow so that the emerging line profile is dominated by the Doppler
    shift caused by the thermal motion of the atoms. The gas temperature defined in the input model
    is used to randomly shift the wavelength of each emitted photon packet according to a
    Maxwell-Boltzmann distribution. This happens outside of this class in the secondary source
    machinery.

    <b>Absorption</b>

    Again following Draine 2011 Chapter 8, and substituting wavelengths for frequencies, the
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
class SpinFlipHydrogenGasMix : public EmittingGasMix
{
    ITEM_CONCRETE(SpinFlipHydrogenGasMix, EmittingGasMix,
                  "A gas mix supporting the spin-flip 21 cm hydrogen transition")
        ATTRIBUTE_TYPE_INSERT(SpinFlipHydrogenGasMix, "CustomMediumState,DynamicState")

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

        PROPERTY_DOUBLE(defaultNeutralSurfaceDensity, "the default neutral hydrogen surface density")
        ATTRIBUTE_QUANTITY(defaultNeutralSurfaceDensity, "masssurfacedensity")
        ATTRIBUTE_MIN_VALUE(defaultNeutralSurfaceDensity, "]0 Msun/pc2")
        ATTRIBUTE_MAX_VALUE(defaultNeutralSurfaceDensity, "1e9 Msun/pc2]")
        ATTRIBUTE_DEFAULT_VALUE(defaultNeutralSurfaceDensity, "10 Msun/pc2")
        ATTRIBUTE_DISPLAYED_IF(defaultNeutralSurfaceDensity, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function determines the radiation field wavelength bin containing the UV field
        strength. */
    void setupSelfBefore() override;

    //======== Capabilities =======

public:
    /** This function returns true, indicating that the cross sections returned by this material
        mix depend on the values of specific state variables other than the number density. */
    bool hasExtraSpecificState() const override;

    /** This function returns DynamicStateType::Secondary, indicating that this material mix has a
        dynamic medium state with updates that affect only secondary emission. */
    DynamicStateType hasDynamicMediumState() const override;

    /** This function returns true, indicating that this material supports secondary line emission
        from gas. */
    bool hasLineEmission() const override;

    //======== Medium state setup =======

public:
    /** This function returns the number and type of import parameters required by this particular
        material mix as a list of SnapshotParameter objects. For this class, the function returns a
        descriptor for the neutral hydrogen mass surface density import parameter. Importing
        metallicity and temperature should be enabled through the corresponding standard
        configuration flags. */
    vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns a list of StateVariable objects describing the specific state
        variables used by the receiving material mix. For this class, the function returns a list
        containing descriptors for the properties defined in the input model (neutral hydrogen
        number density, metallicity, temperature, and neutral hydrogen mass surface density) and
        for a variable to hold the atomic hydrogen fraction derived from the radiation field when
        the dynamic medium state is updated. */
    vector<StateVariable> specificStateVariableInfo() const override;

    /** This function initializes the specific state variables requested by this material mix
        through the specificStateVariableInfo() function except for the number density. For this
        class, the function initializes the temperature, metallicity and neutral hydrogen mass
        surface density to the specified imported values, or if not available, to the
        user-configured default values. The atomic hydrogen fraction is set to zero. */
    void initializeSpecificState(MaterialState* state, double metallicity, double temperature,
                                 const Array& params) const override;

    //======== Medium state updates =======

    /** Based on the specified radiation field and the input model properties found in the given
        material state, this function performs the partitioning scheme described in the class
        header and stores the resulting atomic hydrogen fraction back in the given material state.
        The function returns the update status as described for the UpdateStatus class.
        */
    UpdateStatus updateSpecificState(MaterialState* state, const Array& Jv) const override;

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

    //======== Secondary emission =======

    /** This function returns a list including a single item: the line center of the 21 cm hydrogen
        spinflip transition. */
    Array lineEmissionCenters() const override;

    /** This function returns a list including a single item: the mass of the particle emitting the
        21 cm line, i.e. the hydrogen atom. */
    Array lineEmissionMasses() const override;

    /** This function returns a list including a single item: the 21 cm line luminosity in the
        spatial cell and medium component represented by the specified material state and the
        receiving material mix when it would be embedded in the specified radiation field. */
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
