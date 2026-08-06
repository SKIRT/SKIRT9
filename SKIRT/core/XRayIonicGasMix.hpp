/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef XRAYIONICGASMIX_HPP
#define XRAYIONICGASMIX_HPP

#include "ArrayTable.hpp"
#include "DipolePhaseFunction.hpp"
#include "MaterialMix.hpp"
#include "PhotonPacket.hpp"

////////////////////////////////////////////////////////////////////

/** The XRayIonicGasMix class describes the material properties of a partially ionised gas in the
    X-ray wavelength range, taking into account the effects of photo-absorption, fluorescence,
    recombination, electron scattering, and optionally　resonant scattering. It supports the line 
    production and transfer for the Lyman-series of hydrogen-like ions and the K-series of helium-like 
    ions.　It is largely an extension of the XRayAtomicGasMix class, supporting all ions instead of 
    just the neutral atoms. To avoid the use of this material mix outside of the regime for which it has
    been designed, all cross sections are forced to zero below 4.3 eV and above 500 keV
    (corresponding approximately to a wavelength range from 2.5 pm to 290 nm).

    The class assumes a gas containing a mixture of (ionized) elements with atomic numbers ranging
    from 1 (hydrogen) to 30 (zinc). Including fully ionized ions, this results in a total of 495
    possible ions. Not all ions have to be used, the user can specify what ions are included with
    the \em ions property. The spatial density distribution of the gas is determined by the product
    of the number density and the abundances of the ions in use. The temperature of the gas can
    also be configured by the user as a constant property. In other words, the abundances and the
    temperature are considered to be spatially constant (for a given medium component), while the
    overall density can obviously vary across space as usual.

    <b>Photo-absorption and fluorescence</b>

    Photo-absorption by an ion is the process where the energy of a photon is used to liberate a
    bound electron from one of the electron shells of the atom. This class supports
    photo-absorption by any of the 495 ions in the gas for any of the electron shells, i.e. up to
    the K, L, M, or N shell depending on the ion.

    Fluorescence (in this context) is the process where an electron from a higher energy level
    "falls" into an empty space created by photo-absorption, emitting a new photon with a different
    energy. For each electron shell and for each possible fluorescence transition towards that
    shell, the \em yield defines the probability that such fluorescence event occurs after an
    electron has been liberated in that shell. This class supports the following fluorescent lines,
    all with energies above \f$4.3 \, \text{eV}\f$: K\f$_{\alpha2}\f$, K\f$_{\alpha1}\f$,
    K\f$_{\beta3}\f$, K\f$_{\beta1}\f$, L\f$_{\beta4}\f$, L\f$_{\beta3}\f$, L\f$_{1,2}\f$,
    L\f$_{1,3}\f$, L\f$_{\eta}\f$, L\f$_{l}\f$, L\f$_{\gamma5}\f$, L\f$_{\beta6}\f$,
    L\f$_{\beta1}\f$, L\f$_{\alpha2}\f$, L\f$_{\alpha1}\f$, M\f$_{1,2}\f$, M\f$_{1,3}\f$,
    M\f$_{2,4}\f$, M\f$_{3,4}\f$, M\f$_{3,5}\f$ M\f$_{2}\f$N\f$_{1}\f$, M\f$_{3}\f$N\f$_{1}\f$
    (these transitions are not all the same as the ones in the XRayAtomicGasMix class!).

    Unlike the XRayAtomicGasMix, the fluorescent lines never have an intrinsic shape.

    Because fluorescence only occurs as the result of a photo-absorption event, this class
    implements fluorescence as a form of scattering (where the wavelength of the photon being
    scattered changes) as opposed to emission. This allows both photo-absorption and fluorescence
    to be treated during primary emission. A possible drawback is that the weaker fluorescence
    lines will be represented by a fairly small number of photon packets.

    <b>Recombination</b>

    For hydrogen-like ions and helium-like ions, photo-absorption cannot be followed by fluorescence 
    since the ion will have no excited electrons left. Instead it can be followed by a radiative 
    recombination event either to the ground state or to an excited level, which will then cascade 
    down to the ground state. This class assumes instantaneous recombination and only implements 
    the Lyman-series and He-like K-series photons produced from cascades from the level \f$i\f$ to 
    the ground state. This will thus only result in line emission and no continuum emission. 
    This allows for this process to be modelled the same way as fluorescence.

    <b>Scattering by bound electrons</b>

    Electrons bound to atoms or free in the gas scatter incoming X-ray photons. Currently this
    class can only model electrons as free using inelastic (Compton) scattering. We completely
    ignore elastic (Rayleigh) scattering and inelastic (bound-Compton) scattering. The user can
    select one of three implementations for electron scattering, these are:

    - \em None: ignore scattering by all electrons.

    - \em Free: use free-electron Compton scattering for all electrons. This assumes all electrons
    (bound or free) scatter using the free-electron Compton scattering cross section.

    - \em FreeWithPolarization: use free-electron Compton scattering for all electrons with support
    for polarization.

    This means that all electrons are treated the same way. A material mix with Fe+0 will thus have
    the same scattering cross section as a material mix with Fe+26.

    <b>Resonant scattering</b>

    Resonant scattering is the process where a photon is absorbed and re-emitted by a bound
    electron in an ion, promoting it from the ground state (\f$n=1\f$) to an excited level (\f$n
    \leq 10\f$) and back again. This class implements these transitions for all hydrogen-like ions
    and helium-like ions up to atomic number 30. For hydrogen-like ions, the 18 prominent electric-dipole
    transition lines fine-structure levels are included for each ion (e.g., Ly\f$\alpha_1\f$, 
    Ly\f$\alpha_2\f$, Ly\f$\beta_1\f$, Ly\f$\beta_2\f$, \f$\ldots\f$, Ly\f$\theta_1\f$, Ly\f$\theta_2\f$, 
    Ly\f$\iota_1\f$, Ly\f$\iota_2\f$). The magnetic-dipole transition lines between the metastable level 
    (2s \f$^2S_{1/2}\f$) and the ground level for \f$Z\geq14\f$, where these liens are non-negligible, 
    are also included. For helium-like ions, the He\f$\alpha\f$ \f$w\f$, \f$x\f$, \f$y\f$ and \f$z\f$ lines
    and several prominent higher K-series lines (e.g., He\f$\beta\f$, He\f$\gamma\f$, \f$\ldots\f$, etc.)
    are included. 

    These transitions may occur either as a direct transition from and to the ground state
    or via a radiative cascade. Direct transitions correspond to coherent scattering events, while
    cascades are incoherent and erase information about the initial photon packet. when a cascade
    occurs, we ignore all intermediate photons produced by radiative decays between excited states, 
    keeping only the final transition to the ground state. This means we can also model the cascade 
    as a single scattering event. The branching data, determining the probability of an (in)coherent transition,
    is obtained from the SPEX database Kaastra et al. (2024).

    <b>Configuring the simulation</b>

    In addition to a medium component configured with the material mix represented by this class,
    simulations will usually include primary sources and possibly a dust medium. There is, however,
    no need to include secondary emission, so the simulation mode can be set to "ExtinctionOnly"
    and there is no need to store the radiation field. Only if the user option \em
    resonantScattering is enabled does the simulation mode need to be set to "LyaExtinctionOnly".
    The resulting continuum spectrum and absorption and emission features can be recorded by a
    single instrument configured with a high-resolution wavelength grid, or separate instruments
    can be configured with wavelength grids to resolve specific features of interest.

    Although there is no secondary emission phase, fluorescent emission and lyman recombination,
    which are modelled as scattering, do contribute to the secondary emission spectrum. To study
    these processes seperately from the background continuum, configure the instrument to record
    flux components seperately and consider the secondary emission component, which includes these
    lines in addition to any flux scattered from them. To ensure that low-intensity lines are
    properly included in this flux, set the advanced property \em minWeightReduction in the \c
    PhotonPacketOptions section to a value of \c 1e10 or so. With the default value of \c 1e4,
    low-intensity line photon packets are killed before having a chance to register in the
    instruments.

    The input model must define the spatial distribution of the number density. The use can then
    define the abundances of the ions in the gas in this material mix. The spatial density is
    simply the product of the number density and the abundances of the ions in use.

    If this material mix is associated with a subclass of ImportedMedium, the spatial density
    distribution is read from an input file. In that case, the ski file attributes \em
    importMetallicity, \em importTemperature, and \em importVariableMixParams must be left at
    'false'. For example, if bulk velocities are also imported for this medium component (i.e. \em
    importVelocity is 'true'), the column order would be \f[ ..., n, v_\mathrm{x}, v_\mathrm{y},
    v_\mathrm{z} \f]

    The relative abundances of the present ions in the gas and the temperature of the gas are
    configured in the ski file as constant properties. In other words, the abundances and the
    temperature are considered to be spatially constant (for a given medium component).

    <b>Photo-absorption cross section</b>

    The total photo-absorption cross section per hydrogen atom for this material mix is obtained by
    accumulating the photo-absorption cross sections for all shells and for all individual
    elements, weighted by element abundance, and convolved with a Gaussian profile reflecting each
    element's thermal velocity. Because the abundances and the temperature are fixed, this
    calculation can be performed during setup and the result stored, discretized on a
    high-resolution wavelength grid for later retrieval.

    Verner and Yakovlev (1995, 1996) provide analytic fits to the photo-absorption cross sections
    \f$\sigma_{ph}(E)\f$ as a function of photon energy \f$E\f$ for all the ions up to Z=30:

    \f[\begin{aligned} \sigma_{ph}(E) &= \begin{cases} 0 & E < E_\mathrm{th} \\ \sigma_0 \, F(y) &
    E_\mathrm{th} \le E < E_\mathrm{max}\\ 0 & E_\mathrm{max} \le E \end{cases}, \\ y &=
    \frac{E}{E_0}\\ F(y) &= \left[(y-1)^2+y_{\rm w}^2 \right]y^{-Q} \left(1+ \sqrt{(y/y_{\rm a})}
    \right )^{-P}, \\ Q&=5.5+l-0.5P, \end{aligned} \f] with \f$E_\mathrm{th}\f$ the tabulated
    ionization threshold energy, \f$E_\mathrm{max}\f$ the constant maximum energy (500 keV) for the
    formula to be valid, \f$\sigma_0\f$, \f$E_0\f$, \f$y_{\rm w}\f$, \f$y_{\rm a}\f$ and \f$P\f$
    five tabulated fitting parameters, and \f$l\f$ the subshell orbital quantum number (\f$l=0, 1,
    2, 3\f$ for s, p, d, f orbitals respectively).

    <b>Fluorescence cross section</b>

    The total fluorescence cross section per hydrogen atom for this material mix is obtained
    similarly, but now including only the K, L and M shell photo-absorption cross sections for each
    ion, multiplied by the appropriate fluorescence yields in addition to the element abundance.
    These yields are obtained from Kaastra & Mewe (1993) for all ions up to Z=30, from neutral down
    to B-like.

    <b>Recombination cross section</b>

    The total recombination cross section per hydrogen-like or helium-like ion for this material mix is obtained
    the exact same was as for fluorescence. The yields are now temperature-dependent however, but
    since the temperature is assumed to be constant this is simply read during the setup. The
    temperature-dependent yields are obtained from Mao & Kaastra (2016).

    <b>Electron scattering - cross section and phase function</b>

    As described above, this class provides a free-electron Compton scattering implementation. Note
    that the cross section must be multiplied by the abundance of the corresponding ion. Since we
    assume all electrons to be free, the cross section for an ion with atomic number \f$Z\f$ is
    simply given by \f$Z\,\sigma_\mathrm{C}\f$, where \f$\sigma_\mathrm{C}\f$ is the
    (wavelength-dependent) Compton cross section for a single free electron. The implementation of
    the scattering events is delegated to the ComptonPhaseFunction class; see there for more
    information on the cross section and phase function for free-electron Compton scattering.

    <b>Electron scattering - photon energy shift</b>

    Compton scattering is inelastic, meaning that the photon transfers a fraction of its energy to
    the electron involved in the interaction. The energy shift \f$E'/E\f$ is given by the Compton
    factor \f$C(\theta, E)\f$ defined above. The implementation is delegated to the
    ComptonPhaseFunction class.

    <b>Resonant scattering - cross section and phase function</b>

    The resonant transitions are broadened both thermally and intrinsically. The resulting
    cross section is described by a Voigt profile, the implementation of this is delegated to the
    \em LyUtils and \em VoigtProfile namespaces. The data for the line energies and Einstein A
    coefficients are obtained from the SPEX database Kaastra et al. (2024).

    Resonant scattering in an electric-dipole transition is described as a linear
    combination of an isotropic and a dipole M\"{u}ller matrix. The relative weights
    of these two components are determined by the total angular momentum of the
    lower level, \f$J\f$, and by　\f$\Delta J = J_\mathrm{upper}-J_\mathrm{lower}\f$ \citep{H1947}.
    With the notation used here, \f$E_1\f$ is the weight of the isotropic component
    and \f$E_2\f$ is the weight of the dipole component.

    <table>
    <caption>Weights of the isotropic (\f$E_1\f$) and dipole (\f$E_2\f$)
    components in the M\"{u}ller matrix for resonant scattering.</caption>
    <tr>
      <th>\f$\Delta J\f$</th>
      <th>\f$E_1\f$</th>
      <th>\f$E_2\f$</th>
    </tr>
    <tr>
      <td>\f$1\f$</td>
      <td>\f$\displaystyle \frac{1}{10}
      \frac{3J(6J+7)}{(J+1)(2J+1)}\f$</td>
      <td>\f$\displaystyle \frac{1}{10}
      \frac{(2J+5)(J+2)}{(J+1)(2J+1)}\f$</td>
    </tr>
    <tr>
      <td>\f$0\f$</td>
      <td>\f$\displaystyle \frac{1}{10}
      \frac{3(2J^2+2J+1)}{J(J+1)}\f$</td>
      <td>\f$\displaystyle \frac{1}{10}
      \frac{(2J-1)(2J+3)}{J(J+1)}\f$</td>
    </tr>
    <tr>
      <td>\f$-1\f$</td>
      <td>\f$\displaystyle \frac{1}{10}
      \frac{3(J+1)(6J-1)}{J(2J+1)}\f$</td>
      <td>\f$\displaystyle \frac{1}{10}
      \frac{(2J-3)(J-1)}{J(2J+1)}\f$</td>
    </tr>
    </table>

    For example, the H-like Ly\f$_{\alpha1}\f$ transition has
    \f$J=1/2\f$ and \f$\Delta J=1\f$, and therefore
    \f$E_1=E_2=1/2\f$; the scattered photon is emitted with equal weights
    for the isotropic and dipole components. The H-like Ly\f$_{\alpha2}\f$
    transition has \f$J=1/2\f$ and \f$\Delta J=0\f$, giving
    \f$E_1=1\f$ and \f$E_2=0\f$; the scattered photon is emitted
    isotropically. For the He-like He\f$\alpha\f$ \f$w\f$ and \f$y\f$
    lines, the lower level has \f$J=0\f$ and the upper level has
    \f$J=1\f$, so that \f$E_1=0\f$ and \f$E_2=1\f$; the scattering is
    described entirely by the dipole component.

    The dipole component always includes polarization. Therefore, the
    \em ElectronScattering property does not affect this treatment in any way.
    The implementation of the dipole phase function can be found in the
    \em DipolePhaseFunction class. If \f$E_2=0\f$, the scattered photon is emitted
    isotropically.

    If the scattering is coherent (no-cascade) the wavelength shift is determined by: \f[ \lambda'
    = \lambda \frac{(1 - \boldsymbol{k}_\mathrm{out} \cdot \boldsymbol{v}_\mathrm{atom} / c)}{(1 -
    \boldsymbol{k}_\mathrm{in} \cdot \boldsymbol{v}_\mathrm{atom} / c)} \f] otherwise (cascade) the
    wavelength shift is only determined by the outgoing direction: \f[ \lambda' = \lambda (1 -
    \boldsymbol{k}_\mathrm{out} \cdot \boldsymbol{v}_\mathrm{atom} / c) \f]

    <b>Performing scattering</b>

    The function performing an actual scattering event randomly selects one of the supported
    scattering channels (i.e. scattering by an electron, fluorescent line emission following a
    photo-absorption event or scattering by a resonant line). The relative probabilities for
    these transitions as a function of incoming photon packet wavelength are also calculated during
    setup. The selected transition determines the scattering mechanism. For electrons, Compton
    scattering is used. For fluorescence, the emission direction is isotropic, and the outgoing
    wavelength is the fluorescence wavelength. For the resonant lines, the direction can be
    either isotropic and unpolarized or following a dipole phase function for both direction and
    polarization. The outcome is determined by the total angular momentum of the lower level \f$J\f$
    and the difference in the total angular momentum between the upper and lower levels \f$\Delta J\f$.

    <b>Thermal dispersion</b>

    The thermal dispersion appropriate for a given interaction depends on the (constant)
    temperature configured for this material mix and on the mass of the interacting atom. For
    electron scattering and fluorescence events, the implementation is straightforward. Once a
    channel has been randomly selected, the magnitude of the interacting atom's thermal velocity is
    taken from a precalculated table and a velocity vector is sampled from the Maxwell
    distribution.

    For photo-absorption, the situation is much more involved. In principle, the full cross section
    curve for each ionization transition must be convolved with a Gaussian kernel of appropriate
    width. In practice, the cross section curves are very smooth except for the step at the
    threshold energy. The effect of the convolution is therefore limited to energies near the
    threshold energy, replacing the infinitely steep step by a sigmoid function. Considering that
    the convolution of a step function with a Gaussion is given by the error function, this class
    uses the following approximation. With \f$E_\mathrm{th}\f$ the threshold energy and
    \f$E_\mathrm{s}\f$ the energy dispersion corresponding to the thermal velocity of the atom
    being ionized, the cross section near \f$E_\mathrm{th}\f$ is replaced by

    \f[ \sigma_{ph}'(E) = \frac{\sigma_{ph}(E_\mathrm{th} +2 E_\mathrm{s})} {2} \, \left[ 1+
    \mathrm{erf} ( \frac{E - E_\mathrm{th}} {E_\mathrm{s}} ) \right] \qquad \mathrm{for} \;
    E_\mathrm{th} -2 E_\mathrm{s} < E < E_\mathrm{th} +2 E_\mathrm{s}. \f]

    The sigmoid is scaled by the value of the actual cross section at the end of the interval,
    achieving good continuity at that point in most cases, including all important transitions.

    The thermal broadening present in the resonant scattering is already built into the Voigt
    profile, its implementation can be found in the \em LyUtils and \em VoigtProfile namespaces.
    */
class XRayIonicGasMix : public MaterialMix
{
    /** The enumeration type indicating the implementation used for scattering by electrons. */
    ENUM_DEF(ElectronScattering, None, Free, FreeWithPolarization)
        ENUM_VAL(ElectronScattering, None, "ignore electron")
        ENUM_VAL(ElectronScattering, Free, "use free-electron Compton scattering for all electrons")
        ENUM_VAL(ElectronScattering, FreeWithPolarization,
                 "use free-electron Compton scattering with support for polarization")
    ENUM_END()

    ITEM_CONCRETE(XRayIonicGasMix, MaterialMix,
                  "An ionised gas mix supporting photo-absorption, fluorescence and resonant Lyman scattering for "
                  "X-ray wavelengths")
        ATTRIBUTE_TYPE_INSERT(XRayIonicGasMix, "GasMix")

        PROPERTY_STRING(ions, "the names of the ions for each element (e.g. H,He+,Fe+25,..)")

        PROPERTY_DOUBLE_LIST(abundances, "the abundances of the ions in the same order as the ions property")

        PROPERTY_DOUBLE(temperature, "the temperature of the gas in K")
        ATTRIBUTE_QUANTITY(temperature, "temperature")
        ATTRIBUTE_MIN_VALUE(temperature, "[3")
        ATTRIBUTE_MAX_VALUE(temperature, "1e9]")
        ATTRIBUTE_DEFAULT_VALUE(temperature, "1e4")
        ATTRIBUTE_DISPLAYED_IF(temperature, "Level2")

        PROPERTY_ENUM(electronScattering, ElectronScattering, "implementation of scattering by electrons")
        ATTRIBUTE_DEFAULT_VALUE(electronScattering, "Free")
        ATTRIBUTE_DISPLAYED_IF(electronScattering, "Level3")

        PROPERTY_BOOL(resonantScattering, "enable Lyman resonant scattering for all hydrogen-like ions")
        ATTRIBUTE_DEFAULT_VALUE(resonantScattering, "false")
        ATTRIBUTE_DISPLAYED_IF(resonantScattering, "Level2")
        ATTRIBUTE_RELEVANT_IF(resonantScattering, "simulationModeLyaExtinctionOnly")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function precalculates relevant cross sections and relative contributions over a
        high-resolution wavelength grid. It also stores presistent data that is used during the
        simulation to perform scattering. */
    void setupSelfBefore() override;

    explicit XRayIonicGasMix(SimulationItem* parent, string ions, vector<double> abundances, double temperature,
                             ElectronScattering electronScattering, bool resonantScattering, bool setup);

    /** The destructor destructs the phase function helpers that were created during setup. */
    ~XRayIonicGasMix();

    //======== Private support functions =======

private:
    /** This function returns the index in the private wavelength grid corresponding to the
        specified wavelength. The parameters for converting a wavelength to the appropriate index
        are stored in data members during setup. */
    int indexForLambda(double lambda) const;

    /** This function returns the precalculated thermal velocity of the atom with the specified
        atomic number \f$Z\f$. */
    double vtherm(int Z) const;

    //============= Capabilities =============

public:
    /** This function returns the fundamental material type represented by this material mix, which
        is MaterialType::Gas. */
    MaterialType materialType() const override;

    /** This function returns true if this material mix supports polarization during scattering
        events. For this class, the function returns true if the \em scatterBoundElectrons property
        has been set to \c FreeWithPolarization or is the \em resonantScattering property is set to
        true. */
    bool hasPolarizedScattering() const override;

    /** This function returns true if scattering for this material mix is resonant. For this class,
        the function returns true if the \em resonantScattering property is set to true. */
    bool hasResonantScattering() const override;

    /** This function returns true, indicating that a scattering interaction for this material mix
        may (and usually does) adjust the wavelength of the interacting photon packet. */
    bool hasScatteringDispersion() const override;

    /** This function returns true, indicating that a scattering interaction for this material mix
        may emulate secondary emission. This is used to implement fluorescence as scattering. */
    bool scatteringEmulatesSecondaryEmission() const override;

    //============= Medium state setup =============

    /** This function returns a list of StateVariable objects describing the specific state
        variables used by the receiving material mix. For this class, the function returns just the
        descriptor for the number density. */
    vector<StateVariable> specificStateVariableInfo() const override;

    //============= Low-level material properties =============

    /** This function returns the mass of a hydrogen atom. */
    double mass() const override;

    /** This function returns the absorption (i.e. \f$\sigma_\mathrm{ext} - \sigma_\mathrm{sca}\f$)
        cross section per hydrogen atom at the given wavelength and using the abundances and
        temperature configured for this material mix. */
    double sectionAbs(double lambda) const override;

    /** This function returns the scattering (i.e. electron scattering, fluorescence, resonant
        scattering) cross section per hydrogen atom at the given wavelength and using the
        abundances and temperature configured for this material mix. */
    double sectionSca(double lambda) const override;

    /** This function returns the extinction cross section per hydrogen atom at the given
        wavelength and using the abundances and temperature configured for this material mix. */
    double sectionExt(double lambda) const override;

    //============= High-level photon life cycle =============

    /** This function returns the absorption (i.e. \f$\sigma_\mathrm{ext} - \sigma_\mathrm{sca}\f$)
        opacity at the given wavelength and material state, using the abundances and temperature
        configured for this material mix. The photon packet properties are not used. */
    double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the scattering (i.e. electron scattering, fluorescence, resonant
        scattering) opacity at the given wavelength and material state, using the abundances and
        temperature configured for this material mix. The photon packet properties are not used. */
    double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the extinction opacity at the given wavelength and material state,
        using the abundances and temperature configured for this material mix. The photon packet
        properties are not used. */
    double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

private:
    /** This private function draws a random scattering channel and atom velocity and stores this
        information in the photon packet's scattering information record, unless a previous
        peel-off stored this already. For fluorescence transitions that support a nonzero line
        width, the function also draws a random wavelength from the line shape. For a resonant
        scattering event it determines whether a cascade occrus and if the phase function is
        isotropic or a dipole. */
    void setScatteringInfoIfNeeded(PhotonPacket* pp, const MaterialState* state, const double lambda) const;

public:
    /** This function calculates the contribution of the medium component associated with this
        material mix to the peel-off photon luminosity, polarization state, and wavelength shift
        for the given wavelength, geometry, material state, and photon properties. The
        contributions to the Stokes vector components are stored in the \em I, \em Q, \em U, \em V
        arguments, which are guaranteed to be initialized to zero by the caller. If there is
        wavelength shift, the new wavelength value replaces the incoming value of the \em lambda
        argument.

        The function first calls the private setScatteringInfoIfNeeded() function to establish a
        random scattering channel and atom velocity for this event. In case the selected channel is
        electron scattering, the peel-off bias weight and wavelength shift are determined by the
        Compton phase function and the selected atom velocity. For fluorescence, the peel-off bias
        weight is trivially one because emission is isotropic and unpolarized, and the outgoing
        wavelength is determined by Doppler-shifting the rest wavelength of the selected
        fluorescence transition for the selected atom velocity. For resonant scattering the
        peel-off bias weight and polarization are determined by its isotropic or dipole phase
        function, and the wavelength is shifted accordingly depending on whether the scattering is
        coherent or through a cascade. */
    bool peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs, Direction bfky,
                           const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function performs a scattering event on the specified photon packet in the spatial
        cell and medium component represented by the specified material state and the receiving
        material mix. It first calls the private setScatteringInfoIfNeeded() function to establish
        a random scattering channel and atom velocity for this event.

        In case the selected channel is electron scattering, the outgoing direction and adjusted
        wavelength are determined by the Compton phase function and the selected atom velocity. For
        fluorescence, emission is isotropic, so the outgoing direction is randomly chosen from the
        isotropic distribution. The outgoing wavelength is determined by Doppler-shifting the rest
        wavelength of the selected fluorescence transition for the selected atom velocity. For
        resonant scattering, the emission direction and polarization are determined by its
        isotropic or dipole phase function. The wavelength is shifted depending on whether the
        scattering is coherent or occurs through a cascade. Coherent scattering depends on both the
        incoming and outgoing directions, whereas a cascade removes the incoming photon information
        and depends only on the outgoing direction. */
    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======== Temperature =======

public:
    /** This function returns an indicative temperature of the material mix when it would be
        embedded in a given radiation field. The implementation in this class ignores the radiation
        field and returns the (spatially constant) temperature configured for this material mix. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //======================== Data Members ========================

public:
    // base class for electron scattering helpers (public because we derive from it in anonymous namespace)
    class ScatteringHelper;

private:
    // all of the data below is calculated in the setupSelfBefore(), but is persistent for use after setup to perform scattering

    struct IonParam
    {
        int Z;  // atomic number
        int N;  // number of electrons
    };
    // Compton scattering -> IonParam
    // Fluorescence (+RR) -> FluorescenceParam
    struct FluorescenceParam
    {
        double vth;     // thermal velocity (m/s)
        double lambda;  // wavelength (m)
        double width;   // width (eV)
    };
    // Resonant scattering -> ResonaceLineParam
    struct ResonantParam
    {
        int ionIndex;
        int lineIndex;          // line index
        double vth;             // thermal velocity (m/s)
        double lambda;          // wavelength (m)
        double a;               // Voigt parameter
        double dipoleFraction;  // Probability for dipole phase function
        Array cumBranchingv;    // normalized cumulative branching
    };

    int _numIons;  // number of ions
    int _numFluo;  // number of fluorescence (+RR) transitions
    int _numLine;  // number of resonant scattering transitions

    // persistent data for scattering
    vector<IonParam> _ionParamv;
    vector<FluorescenceParam> _fluorescenceParamv;
    vector<ResonantParam> _resonantParamv;
    vector<double> _vthermv;  // indexed on Z

    // wavelength grid (shifted to the left of the actually sampled points to approximate rounding)
    Array _lambdav;
    // cross sections
    Array _sigmaextv;              // indexed on lambdav
    Array _sigmascav;              // indexed on lambdav
    ArrayTable<2> _cumsigmascavv;  // indexed on lambdav, interaction (electron + fluorescence + resonant)

    // compton-electron scattering helpers depending on the configured implementation
    ScatteringHelper* _com{nullptr};  // Compton scattering helper
    // dipole phase function for resonant scattering
    DipolePhaseFunction* _dpf{nullptr};
};

#endif
