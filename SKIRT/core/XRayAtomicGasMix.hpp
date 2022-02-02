/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef XRAYATOMICGASMIX_HPP
#define XRAYATOMICGASMIX_HPP

#include "ArrayTable.hpp"
#include "MaterialMix.hpp"

////////////////////////////////////////////////////////////////////

/** The XRayAtomicGasMix class describes the material properties related to photo-absorption and
    fluorescence by atomic gas in the X-ray wavelength range. The class assumes a gas containing a
    mixture of non-ionized elements with atomic numbers from 1 (hydrogen) up to 30 (zinc). The
    spatial density distribution of the gas is established by setting the hydrogen density. The
    relative abundances of the 30 elements and the temperature of the gas can be configured by the
    user as constant properties. In other words, the abundances and the temperature are considered
    to be spatially constant (for a given medium component), while the overall density can
    obviously vary across space as usual.

    Photo-absorption by an atom is the process where the energy of a photon is used to liberate a
    bound electron from one of the electron shells of the atom. This class supports
    photo-absorption by any of the 30 elements in the gas for any of the electron shells, i.e. up
    to the K, L, M, or N shell depending on the atomic number of the element. The corresponding
    cross sections can be found in the literature (see below).

    Fluorescence (in this context) is the process where an electron from a higher energy level
    "falls" into an empty space created by photo-absorption, emitting a new photon with a different
    energy. For each electron shell and for each possible fluorescence transition towards that
    shell, the \em yield defines the the probability that such fluorescence event occurs after an
    electron has been liberated in that shell. This class supports K\f$\alpha\f$ and K\f$\beta\f$
    fluorescence (transitions from higher shells towards the K shell) for all elements in the gas.
    The corresponding yields can be found in the literature (see below).

    Because fluorescence only occurs as the result of a photo-absorption event, this class
    implements fluorescence as a form of scattering (where the wavelength of the photon being
    scattered changes) as opposed to emission. This allows both photo-absorption and fluorescence
    to be treated during primary emission. A possible drawback is that the weaker fluorescence
    lines will be represented by a fairly small number of photon packets.

    <b>Configuring the simulation</b>

    In addition to a medium component configured with the material mix represented by this class,
    simulations will usually include primary sources and possibly a dust medium. As indicated
    above, there is no need to include secondary emission, so the simulation mode can be set to
    "ExtinctionOnly" and there is no need to store the radiation field. The resulting continuum
    spectrum and absorption and emission features can be recorded by a single instrument configured
    with a high-resolution wavelength grid, or separate instruments can be configured with
    wavelength grids to resolve specific features of interest.

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

    <b>Extinction (photo-absorption) cross section</b>

    The total extinction cross section per hydrogen atom for this material mix is obtained by
    accumulating the photo-absorption cross sections for all shells and for all individual
    elements, weighted by element abundancy, and convolved with a Gaussian profile reflecting each
    element's thermal velocity. Because the abundancies and the temperature are fixed, this
    calculation can be performed during setup and the result stored, discretized on a
    high-resolution wavelength grid for later retrieval.

    Verner and Yakovlev (1995, www.pa.uky.edu/~verner/photo.html) provide analytic fits to the
    photo-absorption cross sections \f$\sigma_{ph}(E)\f$ as a function of photon energy \f$E\f$ for
    the ground-state shells of the first 30 atomic elements:

    \f[\begin{aligned} \sigma_{ph}(E) &= \begin{cases} 0 & E < E_\mathrm{th} \\ \sigma_0 \,
    F(E/E_0) & E \ge E_\mathrm{th} \end{cases}, \\ F(y) &= \left[(y-1)^2+y_{\rm w}^2 \right]y^{-Q}
    \left(1+ \sqrt{(y/y_{\rm a})} \right )^{-P}, \\ Q&=5.5+l-0.5P, \end{aligned} \f]

    with \f$E_\mathrm{th}\f$ the tabulated ionization threshold energy, \f$\sigma_0\f$, \f$E_0\f$,
    \f$y_{\rm w}\f$, \f$y_{\rm a}\f$ and \f$P\f$ five tabulated fitting parameters, and \f$l\f$ the
    subshell orbital quantum number (\f$l=0, 1, 2, 3\f$ for s, p, d, f orbitals respectively).

    <b>Scattering (fluorescence) cross section</b>

    The total "scattering" cross section per hydrogen atom for this material mix is obtained
    similarly, but now including only the K shell photo-absorption cross section for each element
    and multiplying by the appropriate fluorescence yields in addition to element abundancy.

    The total "absorption" cross section is then simply obtained by subtracting the "scattering"
    cross section from the extinction cross section.

    <b>Performing scattering (fluorescence)</b>

    The function performing an actual scattering event randomly selects one of the supported
    fluorescence transitions (i.e. K\f$\alpha\f$ or K\f$\beta\f$ for one of the supported
    elements). The relative probabilities for these transitions (as a function of incoming photon
    packet wavelength) are also calculated during setup. The selected transition determines the
    fluorescence wavelength. The outgoing photon packet wavelength is then obtained by adding a
    random Gaussian dispersion reflecting the interacting element's thermal velocity. The emission
    direction is isotropic.

    */
class XRayAtomicGasMix : public MaterialMix
{
    ITEM_CONCRETE(XRayAtomicGasMix, MaterialMix,
                  "A gas mix supporting photo-absorption and fluorescence for X-ray wavelengths")

        PROPERTY_DOUBLE_LIST(abundancies, "the abundancies for the elements with atomic number Z = 1,...,30")
        ATTRIBUTE_MIN_VALUE(abundancies, "[0")
        ATTRIBUTE_MAX_VALUE(abundancies, "1]")
        ATTRIBUTE_REQUIRED_IF(abundancies, "false")
        ATTRIBUTE_DISPLAYED_IF(abundancies, "Level2")

        PROPERTY_DOUBLE(temperature, "the temperature of the gas")
        ATTRIBUTE_QUANTITY(temperature, "temperature")
        ATTRIBUTE_MIN_VALUE(temperature, "[3")  // gas temperature must be above local Universe T_CMB
        ATTRIBUTE_MAX_VALUE(temperature, "1e9]")
        ATTRIBUTE_DEFAULT_VALUE(temperature, "1e4")
        ATTRIBUTE_DISPLAYED_IF(temperature, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function precalculates relevant cross sections and relative constribution over a
        high-resolution wavelength grid. */
    void setupSelfBefore() override;

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
        may adjust the wavelength of the interacting photon packet. */
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

    /** This function returns the absorption (i.e. extinction minus fluorescence) cross section per
        hydrogen atom at the given wavelength and using the abundances and temperature configured
        for this material mix. */
    double sectionAbs(double lambda) const override;

    /** This function returns the scattering (i.e. fluorescence) cross section per hydrogen atom at
        the given wavelength and using the abundances and temperature configured for this material
        mix. */
    double sectionSca(double lambda) const override;

    /** This function returns the extinction cross section per hydrogen atom at the given
        wavelength and using the abundances and temperature configured for this material mix. */
    double sectionExt(double lambda) const override;

    //======== High-level photon life cycle =======

    /** This function returns the absorption (i.e. extinction minus fluorescence) opacity at the
        given wavelength and material state, using the abundances and temperature configured for
        this material mix. The photon packet properties are not used. */
    double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the scattering (i.e. fluorescence) opacity at the given wavelength
        and material state, using the abundances and temperature configured for this material mix.
        The photon packet properties are not used. */
    double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the extinction opacity at the given wavelength and material state,
        using the abundances and temperature configured for this material mix. The photon packet
        properties are not used. */
    double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function calculates the contribution of the medium component associated with this
        material mix to the peel-off photon luminosity, polarization state, and wavelength shift
        for the given wavelength, geometry, material state, and photon properties. The
        contributions to the Stokes vector components are added to the incoming values of the \em
        I, \em Q, \em U, \em V arguments. If there is wavelength shift, the new wavelength value
        replaces the incoming value of the \em lambda argument.

        Since we force the peel-off photon packet to be scattered from the direction \f${\bf{k}}\f$
        into the direction \f${\bf{k}}_{\text{obs}}\f$, the corresponding biasing factor is given
        by the probability that a photon packet would be scattered into the direction
        \f${\bf{k}}_{\text{obs}}\f$ if its original propagation direction was \f${\bf{k}}\f$. For a
        given medium component, this biasing factor is equal to the value of the scattering phase
        function \f$\Phi({\bf{k}},{\bf{k}}_{\text{obs}})\f$ for that medium component. If there are
        multiple medium components, the aggregated biasing factor is the mean of the scattering
        phase function values weighted using the relative opacities for the various components. The
        relative opacity weight for the current component is specified as argument \em w. */
    void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, double w, Direction bfkobs,
                           Direction bfky, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function performs a scattering event on the specified photon packet in the spatial
        cell and medium component represented by the specified material state and the receiving
        material mix. Most of the properties of the photon packet remain unaltered, including the
        position and the luminosity. The properties that change include the number of scattering
        events experienced by the photon packet, which is increased by one, the propagation
        direction, which is generated randomly, the wavelength, which is properly Doppler-shifted
        for the bulk velocity of the medium, and the polarization state, which may be affected by
        the scattering process.

        The calculation takes all physical processes into account, including the bulk velocity and
        Hubble expansion velocity in the cell, any relevant material state variables such as the
        temperature of a gas medium, and any relevant properties of the incoming photon packet such
        as the polarization state. The first argument specifies the perceived wavelength of the
        photon packet at the scattering location so that this value does not need to be
        recalculated within the function. */
    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======== Temperature =======

    /** This function returns an indicative temperature of the material mix when it would be
        embedded in a given radiation field. The implementation in this class ignores the radiation
        field and returns the (spatially constant) temperature configured for this material mix. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //======================== Data Members ========================

private:
    // all data members are precalculated in setupSelfAfter()

    // wavelength grid (shifted to the left of the actually sampled points to approximate rounding)
    Array _lambdav;  // indexed on ell

    // cross sections
    Array _sigmaextv;  // indexed on ell
    Array _sigmascav;  // indexed on ell

    // emission wavelengths and normalized cumulative probability distributions of fluorescence channels
    vector<double> _fluolambdav;   // indexed on k
    ArrayTable<2> _fluocumprobvv;  // indexed on ell, k
};

////////////////////////////////////////////////////////////////////

#endif
