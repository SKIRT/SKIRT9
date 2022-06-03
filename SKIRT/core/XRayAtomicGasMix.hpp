/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef XRAYATOMICGASMIX_HPP
#define XRAYATOMICGASMIX_HPP

#include "ArrayTable.hpp"
#include "MaterialMix.hpp"
#include "PhotonPacket.hpp"

////////////////////////////////////////////////////////////////////

/** The XRayAtomicGasMix class describes the material properties of neutral atomic gas in the X-ray
    wavelength range, taking into account the effects of photo-absorption, fluorescence, and
    scattering by bound electrons. To avoid the use of this material mix outside of the regime for
    which it has been designed, all cross sections are forced to zero below 4.3 eV and above 500
    keV (corresponding approximately to a wavelength range from 2.5 pm to 290 nm).

    The class assumes a gas containing a mixture of non-ionized elements with atomic numbers from 1
    (hydrogen) up to 30 (zinc). The spatial density distribution of the gas is established by
    setting the hydrogen density. The relative abundances of the 30 elements and the temperature of
    the gas can be configured by the user as constant properties. In other words, the abundances
    and the temperature are considered to be spatially constant (for a given medium component),
    while the overall density can obviously vary across space as usual.

    Photo-absorption by an atom is the process where the energy of a photon is used to liberate a
    bound electron from one of the electron shells of the atom. This class supports
    photo-absorption by any of the 30 elements in the gas for any of the electron shells, i.e. up
    to the K, L, M, or N shell depending on the atomic number of the element.

    Fluorescence (in this context) is the process where an electron from a higher energy level
    "falls" into an empty space created by photo-absorption, emitting a new photon with a different
    energy. For each electron shell and for each possible fluorescence transition towards that
    shell, the \em yield defines the probability that such fluorescence event occurs after an
    electron has been liberated in that shell. This class supports K\f$\alpha\f$ and K\f$\beta\f$
    fluorescence (transitions from higher shells towards the K shell) for all elements in the gas.

    Because fluorescence only occurs as the result of a photo-absorption event, this class
    implements fluorescence as a form of scattering (where the wavelength of the photon being
    scattered changes) as opposed to emission. This allows both photo-absorption and fluorescence
    to be treated during primary emission. A possible drawback is that the weaker fluorescence
    lines will be represented by a fairly small number of photon packets.

    Electrons bound to the atoms in the gas scatter incoming X-ray photons. This process can be
    elastic (Rayleigh scattering) or inelastic (bound-Compton scattering). The user can select one
    of four implementations of bound-electron scattering. In increasing order of accuracy and
    computational effort, these are:

    - \em None: ignore scattering by bound electrons. This approximation works reasonably well for
    low energies only (as the total extinction cross section is dominated by photo-absorption).
    However, it is wrong for higher energies. This models photo-absorption and fluorescence only.

    - \em Free: use free-electron Compton scattering. This approximation works reasonably well for
    high energies (because the bound-electron scattering cross section approaches that of
    free-electron scattering) and for lower energies (because the total scattering cross section is
    dominated by photo-absorption anyway). However, it is less accurate for intermediate energies.

    - \em Good: use the smooth Rayleigh scattering approximation and exact bound-Compton
    scattering. This reflects the implementation by most X-ray codes and can be considered to be a
    good approximation. It relies on tabulated cross sections and correction factors.

    - \em Exact: use anomalous Rayleigh scattering and exact bound-Compton scattering. This
    implementation accurately reflects the physical processes, apart from approximations caused by
    discretization, tabulation, and interpolation.

    <b>Configuring the simulation</b>

    In addition to a medium component configured with the material mix represented by this class,
    simulations will usually include primary sources and possibly a dust medium. As indicated
    above, there is no need to include secondary emission, so the simulation mode can be set to
    "ExtinctionOnly" and there is no need to store the radiation field. The resulting continuum
    spectrum and absorption and emission features can be recorded by a single instrument configured
    with a high-resolution wavelength grid, or separate instruments can be configured with
    wavelength grids to resolve specific features of interest.

    To study fluorescence lines seperated from the background continuum, configure the instrument
    to record flux components seperately and consider the scattered flux component, which includes
    the fluorescence lines in addition to any flux scattered from other media components. To ensure
    that low-intensity lines are properly included in this flux, set the advanced property \em
    minWeightReduction in the \c PhotonPacketOptions section to a value of \c 1e10 or so. With the
    default value of \c 1e4, low-intensity fluorescence photon packets are killed before having a
    chance to register in the instruments.

    The input model must define the spatial distribution of the hydrogen number density
    \f$n_\mathrm{H} = n_\mathrm{HI} + n_\mathrm{HII} + 2\,n_\mathrm{H2}\f$, i.e. including atomic,
    ionized, and molecular hydrogen. Note the factor 2 in this equation; \f$n_\mathrm{H}\f$ in fact
    specifies the number of hydrogen protons per volume rather than the number of hydrogen-like
    particles per volume.

    If this material mix is associated with a geometric medium component, the geometry defines the
    spatial density distribution \f$n_\mathrm{H}\f$. Alternatively, if the material mix is
    associated with a subclass of ImportedMedium, the spatial density distribution is read from an
    input file. In that case, the ski file attributes \em importMetallicity, \em importTemperature,
    and \em importVariableMixParams must be left at 'false'. For example, if bulk velocities are
    also imported for this medium component (i.e. \em importVelocity is 'true'), the column order
    would be \f[ ..., n_\mathrm{H}, v_\mathrm{x}, v_\mathrm{y}, v_\mathrm{z} \f]

    The relative abundances of the 30 elements in the gas and the temperature of the gas are
    configured in the ski file as constant properties. In other words, the abundances and the
    temperature are considered to be spatially constant (for a given medium component). If the list
    of abundances is left empty in the ski file, default abundancies are used, taken from Table 2
    of Anders & Grevesse (1989), the default abundance table in Xspec. In the default list, the
    abundance of hydrogen is set to unity. However, it is acceptable to specify a hydrogen
    abundance lower than one, for example to model an ionized hydrogen fraction.

    <b>Photo-absorption cross section</b>

    The total photo-absorption cross section per hydrogen atom for this material mix is obtained by
    accumulating the photo-absorption cross sections for all shells and for all individual
    elements, weighted by element abundancy, and convolved with a Gaussian profile reflecting each
    element's thermal velocity. Because the abundancies and the temperature are fixed, this
    calculation can be performed during setup and the result stored, discretized on a
    high-resolution wavelength grid for later retrieval.

    Verner and Yakovlev (1995, 1996) provide analytic fits to the photo-absorption cross sections
    \f$\sigma_{ph}(E)\f$ as a function of photon energy \f$E\f$ for the ground-state shells of the
    first 30 atomic elements:

    \f[\begin{aligned} \sigma_{ph}(E) &= \begin{cases} 0 & E < E_\mathrm{th} \\ \sigma_0 \,
    F(y) & E_\mathrm{th} \le E < E_\mathrm{max}\\ 0 & E_\mathrm{max} \le E \end{cases}, \\
    x &= \frac{E}{E_0} - y_0, \; \; \; \; y = \sqrt{x^2 + y_1^2},\\ F(y) &= \left[(x-1)^2+y_{\rm w}^2
    \right]y^{-Q} \left(1+ \sqrt{(y/y_{\rm a})} \right )^{-P}, \\ Q&=5.5+l-0.5P, \end{aligned} \f]
    with \f$E_\mathrm{th}\f$ the tabulated ionization threshold energy, \f$E_\mathrm{max}\f$ the
    tabulated maximum energy for the formula to be valid, \f$\sigma_0\f$, \f$E_0\f$, \f$y_{\rm w}\f$,
    \f$y_{\rm a}\f$, \f$P\f$, \f$y_0\f$ and \f$y_1\f$ seven tabulated fitting parameters, and \f$l\f$
    the subshell orbital quantum number (\f$l=0, 1, 2, 3\f$ for s, p, d, f orbitals respectively).

    <b>Fluorescence cross section</b>

    The total fluorescence cross section per hydrogen atom for this material mix is obtained
    similarly, but now including only the K shell photo-absorption cross section for each element
    and multiplying by the appropriate fluorescence yields in addition to element abundancy.

    <b>Electron scattering</b>

    As described above, this class provides four implementations for the scattering of X-rays
    photons by the electrons bound to the atoms. These implementations involve four types of
    scattering: free-electron Compton scattering, bound-electron Compton scattering, smooth
    Rayleigh scattering, and anomalous Rayleigh scattering. Note that, in all cases, the listed
    cross sections must be multiplied by the abundance of the corresponding element relative to
    hydrogen.

    For free-electron Compton scattering, the cross section for an element with atomic number
    \f$Z\f$ is given by \f$Z\,\sigma_\mathrm{C}\f$, where \f$\sigma_\mathrm{C}\f$ is the
    (wavelength-dependent) Compton cross section for a single free electron. The implementation of
    the scattering events is delegated to the ComptonPhaseFunction class; see there for more
    information on the cross section and phase function for free-electron Compton scattering.

    The other scattering implementations depend on quantities that are available for each element
    \f$Z\f$ in tabulated form as a function of the incoming photon energy \f$E\f$ or as a function
    of the dimensionless momentum transfer parameter \f$q\f$, \f[q= \frac{E}{12.4 \, \mathrm{keV}}
    \cdot \sin(\theta/2), \f] with \f$\theta\f$ the scattering angle.

    For bound-electron Compton scattering, the cross sections \f$\sigma_{CS, Z}(E)\f$ are available
    as a table. The normalised scattering phase function for element Z is given by \f[ \Phi_{CS,
    Z}(\theta, E)= \frac{3}{4}\, \frac{\sigma_T}{\sigma_{CS, Z}(E)}\Big[C^3(\theta, E) +
    C(\theta, E) -C^2(\theta, E)\sin^2\theta\Big] \cdot S_Z(q), \f] with tabulated incoherent
    scattering functions \f$S_Z(q)\f$ and the Compton factor \f$C(\theta, E)\f$ defined as \f[
    C(\theta, E) = {\Big[{1+\frac{E}{m_ec^2}(1-\cos \theta)\Big]}}^{-1}. \f] Also, inelastic
    bound-Compton scattering will change the photon energy by the Compton factor, just as for
    free-electron scattering.

    For smooth Rayleigh scattering, the cross sections \f$\sigma_{RSS, Z}(E)\f$ are available as a
    table. The normalised scattering phase function for element Z is given by \f[ \Phi_{RSS,
    Z}(\theta, E)= \frac{3}{4}\, \frac{\sigma_T}{\sigma_{RSS, Z}(E)}\Big[ 1 + \cos^2\theta
    \Big] \cdot F_Z^2(q), \f] with tabulated atomic form factors \f$F_Z(q)\f$, which converge to
    \f$Z\f$ at small \f$q\f$ and decrease to zero for large \f$q\f$.

    Similarly, for anomalous Rayleigh scattering, the cross sections \f$\sigma_{RSA, Z}(E)\f$ are
    available as a table. The normalised scattering phase function for element Z is now given by
    \f[ \Phi_{RSA, Z}(\theta, E)= \frac{3}{4}\, \frac{\sigma_T}{\sigma_{RSA, Z}(E)}\Big[ 1 +
    \cos^2\theta \Big] \cdot \Big[\big(F_Z(q) + F'_Z(E)\big)^2 + {F''_Z}^2(E)\Big], \f] with the
    same atomic form factors \f$F_Z(q)\f$ as before, and tabulated real and imaginary anomalous
    scattering functions \f$F'_Z(E)\f$ and \f$F''_Z(E)\f$.

    <b>Performing scattering</b>

    The function performing an actual scattering event randomly selects one of the supported
    scattering channels (i.e. scattering by an electron bound to one of the supported elements or a
    K\f$\alpha\f$ or K\f$\beta\f$ fluorescence transition for one of the supported elements). The
    relative probabilities for these transitions as a function of incoming photon packet wavelength
    are also calculated during setup. The selected transition determines the scattering mechanism.
    For bound electrons, Rayleigh or Compton scattering is used. For fluorescence, the emission
    direction is isotropic, and the outgoing wavelength is the fluorescence wavelength. In both
    cases, a random Gaussian dispersion reflecting the interacting element's thermal velocity is
    applied to the outgoing wavelength.

    <b>Thermal dispersion</b>

    The thermal dispersion appropriate for a given interaction depends on the (constant)
    temperature configured for this material mix and on the mass of the interacting atom. For
    electron scattering and fluorescence events, the implementation is straightforward. Once a
    channel has been randomly selected, the magnitude of the interacting atom's thermal velocity is
    taken from a precalculated table and a velocity vector is sampled from the Maxwell
    distribution.

    For absorption, the situation is much more involved. In principle, the full cross section curve
    for each ionization transition must be convolved with a Gaussian kernel of appropriate width.
    In practice, the cross section curves are very smooth except for the step at the threshold
    energy. The effect of the convolution is therefore limited to energies near the threshold
    energy, replacing the infinitely steep step by a sigmoid function. Considering that the
    convolution of a step function with a Gaussion is given by the error function, this class uses
    the following approximation. With \f$E_\mathrm{th}\f$ the threshold energy and
    \f$E_\mathrm{s}\f$ the energy dispersion corresponding to the thermal velocity of the atom
    being ionized, the cross section near \f$E_\mathrm{th}\f$ is replaced by

    \f[ \sigma_{ph}'(E) = \frac{\sigma_{ph}(E_\mathrm{th} +2 E_\mathrm{s})} {2} \, \left[ 1+
    \mathrm{erf} ( \frac{E - E_\mathrm{th}} {E_\mathrm{s}} ) \right] \qquad \mathrm{for} \;
    E_\mathrm{th} -2 E_\mathrm{s} < E < E_\mathrm{th} +2 E_\mathrm{s}. \f]

    The sigmoid is scaled by the value of the actual cross section at the end of the interval,
    achieving good continuity at that point in most cases, including all important transitions.

    */
class XRayAtomicGasMix : public MaterialMix
{
    /** The enumeration type indicating the implementation used for scattering by bound electrons.
        */
    ENUM_DEF(BoundElectrons, None, Free, Good, Exact)
        ENUM_VAL(BoundElectrons, None, "ignore bound electrons")
        ENUM_VAL(BoundElectrons, Free, "use free-electron Compton scattering")
        ENUM_VAL(BoundElectrons, Good, "use smooth Rayleigh scattering and exact bound-Compton scattering")
        ENUM_VAL(BoundElectrons, Exact, "use anomalous Rayleigh scattering and exact bound-Compton scattering")
    ENUM_END()

    ITEM_CONCRETE(XRayAtomicGasMix, MaterialMix,
                  "A gas mix supporting photo-absorption and fluorescence for X-ray wavelengths")
        ATTRIBUTE_TYPE_INSERT(XRayAtomicGasMix, "GasMix")

        PROPERTY_DOUBLE_LIST(abundancies, "the abundancies for the elements with atomic number Z = 1,...,30")
        ATTRIBUTE_MIN_VALUE(abundancies, "[0")
        ATTRIBUTE_MAX_VALUE(abundancies, "1]")
        ATTRIBUTE_REQUIRED_IF(abundancies, "false")
        ATTRIBUTE_DISPLAYED_IF(abundancies, "Level2")

        PROPERTY_DOUBLE(temperature, "the temperature of the gas or zero to disable thermal dispersion")
        ATTRIBUTE_QUANTITY(temperature, "temperature")
        ATTRIBUTE_MIN_VALUE(temperature, "[0")
        ATTRIBUTE_MAX_VALUE(temperature, "1e9]")
        ATTRIBUTE_DEFAULT_VALUE(temperature, "1e4")
        ATTRIBUTE_DISPLAYED_IF(temperature, "Level2")

        PROPERTY_ENUM(scatterBoundElectrons, BoundElectrons, "implementation of scattering by bound electrons")
        ATTRIBUTE_DEFAULT_VALUE(scatterBoundElectrons, "Good")
        ATTRIBUTE_DISPLAYED_IF(scatterBoundElectrons, "Level3")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function precalculates relevant cross sections and relative contributions over a
        high-resolution wavelength grid. */
    void setupSelfBefore() override;

    /** The destructor destructs the phase function helpers that were created during setup. */
    ~XRayAtomicGasMix();

    //======== Private support functions =======

private:
    /** This function returns the index in the private wavelength grid corresponding to the
        specified wavelength. The parameters for converting a wavelength to the appropriate index
        are stored in data members during setup. */
    int indexForLambda(double lambda) const;

    //======== Capabilities =======

public:
    /** This function returns the fundamental material type represented by this material mix, which
        is MaterialType::Gas. */
    MaterialType materialType() const override;

    /** This function returns true, indicating that a scattering interaction for this material mix
        may (and usually does) adjust the wavelength of the interacting photon packet. */
    bool hasScatteringDispersion() const override;

    //======== Medium state setup =======

public:
    /** This function returns a list of StateVariable objects describing the specific state
        variables used by the receiving material mix. For this class, the function returns just the
        descriptor for the number density. */
    vector<StateVariable> specificStateVariableInfo() const override;

    //======== Low-level material properties =======

public:
    /** This function returns the mass of a hydrogen atom. */
    double mass() const override;

    /** This function returns the absorption (i.e. extinction minus scattering) cross section per
        hydrogen atom at the given wavelength and using the abundances and temperature configured
        for this material mix. */
    double sectionAbs(double lambda) const override;

    /** This function returns the scattering (i.e. bound electrons plus fluorescence) cross section
        per hydrogen atom at the given wavelength and using the abundances and temperature
        configured for this material mix. */
    double sectionSca(double lambda) const override;

    /** This function returns the extinction cross section per hydrogen atom at the given
        wavelength and using the abundances and temperature configured for this material mix. */
    double sectionExt(double lambda) const override;

    //======== High-level photon life cycle =======

public:
    /** This function returns the absorption (i.e. extinction minus scattering) opacity at the
        given wavelength and material state, using the abundances and temperature configured for
        this material mix. The photon packet properties are not used. */
    double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the scattering (i.e. bound electrons plus fluorescence) opacity at
        the given wavelength and material state, using the abundances and temperature configured
        for this material mix. The photon packet properties are not used. */
    double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the extinction opacity at the given wavelength and material state,
        using the abundances and temperature configured for this material mix. The photon packet
        properties are not used. */
    double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

private:
    /** This private function draws a random scattering channel and atom velocity and stores this
        information in the photon packet's scattering information record, unless a previous
        peel-off stored this already. */
    void setScatteringInfoIfNeeded(PhotonPacket::ScatteringInfo* scatinfo, double lambda) const;

public:
    /** This function calculates the contribution of the medium component associated with this
        material mix to the peel-off photon luminosity and wavelength shift for the given
        wavelength and observer direction. It first calls the private setScatteringInfoIfNeeded()
        function to establish a random scattering channel and atom velocity for this event.

        In case the selected channel is bound electron scattering, the peel-off bias weight and
        wavelength shift are determined by the Compton phase function and the selected atom
        velocity. For fluorescence, the peel-off bias weight is trivially one because emission is
        isotropic, and the outgoing wavelength is determined by Doppler-shifting the rest
        wavelength of the selected fluorescence transition for the selected atom velocity. */
    void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs, Direction bfky,
                           const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function performs a scattering event on the specified photon packet in the spatial
        cell and medium component represented by the specified material state and the receiving
        material mix. It first calls the private setScatteringInfoIfNeeded() function to establish
        a random scattering channel and atom velocity for this event.

        In case the selected channel is bound electron scattering, the outgoing direction and
        adjusted wavelength are determined by the Compton phase function and the selected atom
        velocity. For fluorescence, emission is isotropic, so the outgoing direction is randomly
        chosen from the isotropic distribution. The outgoing wavelength is determined by
        Doppler-shifting the rest wavelength of the selected fluorescence transition for the
        selected atom velocity. */
    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======== Temperature =======

public:
    /** This function returns an indicative temperature of the material mix when it would be
        embedded in a given radiation field. The implementation in this class ignores the radiation
        field and returns the (spatially constant) temperature configured for this material mix. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //======================== Data Members ========================

public:
    // base class for bound-electron scattering helpers (public because we derive from it in anonymous namespace)
    class ScatteringHelper;

private:
    // all data members are precalculated in setupSelfAfter()

    // wavelength grid (shifted to the left of the actually sampled points to approximate rounding)
    Array _lambdav;  // indexed on ell

    // total extinction and scattering cross sections
    Array _sigmaextv;  // indexed on ell
    Array _sigmascav;  // indexed on ell

    // emission wavelengths for each of the fluorescence transitions
    vector<double> _lambdafluov;  // indexed on k

    // thermal velocities and normalized cumulative probability distributions for the scattering channnels:
    //   - Rayleigh scattering by bound electrons for each atom
    //   - Compton scattering by bound electrons for each atom
    //   - fluorescence transitions
    vector<double> _vthermscav;   // indexed on m
    ArrayTable<2> _cumprobscavv;  // indexed on ell, m

    // bound-electron scattering helpers depending on the configured implementation
    ScatteringHelper* _ray{nullptr};  // Rayleigh scattering helper
    ScatteringHelper* _com{nullptr};  // Compton scattering helper
};

////////////////////////////////////////////////////////////////////

#endif
