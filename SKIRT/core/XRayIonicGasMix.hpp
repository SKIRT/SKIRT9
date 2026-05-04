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
    Lyman recombination, electron scattering, and optionally resonant Lyman resonant scattering. It
    is largely an extension of the XRayAtomicGasMix class, supporting all ions instead of just the
    neutral atoms. To avoid the use of this material mix outside of the regime for which it has
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

    <b>Lyman recombination continuum</b>

    For hydrogen-like ions, photo-absorption can not be followed by fluorescence since the ion will
    have no electrons left. Instead it can be followed by a radiative recombination event either to
    the ground state or to an excited level, which will then cascade down to the ground state. This
    class only implements the Lyman series photons as a product of the cascade from the level
    \f$i\f$ to the ground state. This will thus only result in line emission and can also be
    modelled the same way as fluorescence.

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

    <b>Performing scattering</b>

    The function performing an actual scattering event randomly selects one of the supported
    scattering channels (i.e. scattering by an electron, fluorescent line emission following
    a photo-absorption event or scattering by a resonant Lyman line). The relative
    probabilities for these transitions as a function of incoming photon packet wavelength are also
    calculated during setup. The selected transition determines the scattering mechanism. For
    electrons, Compton scattering is used. For fluorescence, the emission direction is
    isotropic, and the outgoing wavelength is the fluorescence wavelength. For the resonant Lyman lines,
    the direction can be either isotropic and unpolarized or following a dipole phase function for 
    both direction and polarization. The outcome is determined by the Lyman index... lya1, lyb1, etc.

    
    


    */
class XRayIonicGasMix : public MaterialMix
{
    ENUM_DEF(ElectronScattering, None, Free, FreeWithPolarization)
        ENUM_VAL(ElectronScattering, None, "ignore electron")
        ENUM_VAL(ElectronScattering, Free, "use free-electron Compton scattering for all electrons")
        ENUM_VAL(ElectronScattering, FreeWithPolarization,
                 "use free-electron Compton scattering with support for polarization")
    ENUM_END()

    ITEM_CONCRETE(XRayIonicGasMix, MaterialMix, "Ionised gas mix")
        ATTRIBUTE_TYPE_INSERT(XRayIonicGasMix, "GasMix,CustomMediumState")

        PROPERTY_STRING(ions, "the names of the ions for each element (e.g. H,He+,Li+1,..)")

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
    explicit XRayIonicGasMix(SimulationItem* parent, string ions, vector<double> abundances, double temperature,
                             ElectronScattering boundElectrons, bool resonantScattering, bool setup);

    void setupSelfBefore() override;

    ~XRayIonicGasMix();

    //======== Private support functions =======

private:
    /** This function returns the index in the private wavelength grid corresponding to the
        specified wavelength. The parameters for converting a wavelength to the appropriate index
        are stored in data members during setup. */
    int indexForLambda(double lambda) const;

    double vtherm(int Z) const;

    //============= Capabilities =============

public:
    MaterialType materialType() const override;

    bool hasPolarizedScattering() const override;

    bool hasResonantScattering() const override;

    bool hasScatteringDispersion() const override;

    bool scatteringEmulatesSecondaryEmission() const override;

    //============= Medium state setup =============

    vector<StateVariable> specificStateVariableInfo() const override;

    //============= Low-level material properties =============

    double mass() const override;

    double sectionAbs(double lambda) const override;

    double sectionSca(double lambda) const override;

    double sectionExt(double lambda) const override;

    //============= High-level photon life cycle =============

    double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    Array sigmaSca(double lambda, const MaterialState* state) const;

    void setScatteringInfoIfNeeded(PhotonPacket* pp, const MaterialState* state, const double lambda) const;

    bool peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs, Direction bfky,
                           const MaterialState* state, const PhotonPacket* pp) const override;

    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======== Temperature =======

public:
    /** This function returns an indicative temperature of the material mix when it would be
        embedded in a given radiation field. The implementation in this class ignores the radiation
        field and returns the (spatially constant) temperature configured for this material mix. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //======================== Data Members ========================

public:
    // base class for bound1-electron scattering helpers (public because we derive from it in anonymous namespace)
    class ScatteringHelper;

private:
    // all data is calculated in the setupSelfBefore(), but is persistent for use after setup to perform scattering

    struct IonParam
    {
        IonParam(int Z, int N) : Z(Z), N(N) {}

        int Z;  // atomic number
        int N;  // number of electrons
    };
    // Compton scattering -> IonParam
    // Fluorescence (+Lyman RC) -> FluorescenceParam
    struct FluorescenceParam
    {
        int Z;          // atomic number
        double lambda;  // wavelength (m)
        double width;   // width (eV)
    };
    // Resonant scattering -> LymanParam
    struct LymanParam
    {
        int Z;                // atomic number
        int index;            // Lyman index (alpha1/2, alpha3/2, beta1/2, ...)
        double lambda;        // wavelength (m)
        double a;             // Voigt parameter
        Array cumbranchingv;  // normalized cumulative branching
    };

    int _numIons;  // number of ions
    int _numFluo;  // number of fluorescence (+Lyman RC) transitions
    int _numLym;   // number of Lyman resonant scattering transitions

    // persistent data for scattering
    vector<IonParam> _ionParamv;
    vector<FluorescenceParam> _fluorescenceParamv;
    vector<LymanParam> _lymanParamv;
    vector<double> _vthermv;  // indexed on Z

    // wavelength grid (shifted to the left of the actually sampled points to approximate rounding)
    Array _lambdav;
    // cross sections
    Array _sigmaextv;              // indexed on lambdav
    Array _sigmascav;              // indexed on lambdav
    ArrayTable<2> _cumsigmascavv;  // indexed on lambdav, interactions (numIons + numFluo + numRes)

    // compton-electron scattering helpers depending on the configured implementation
    ScatteringHelper* _com{nullptr};  // Compton scattering helper
    // dipole phase function for resonant scattering
    DipolePhaseFunction* _dpf{nullptr};
};

#endif
