/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTMIX_HPP
#define DUSTMIX_HPP

#include "ArrayTable.hpp"
#include "EquilibriumDustEmissionCalculator.hpp"
#include "MaterialMix.hpp"
#include "Range.hpp"
#include "Table.hpp"

////////////////////////////////////////////////////////////////////

/** DustMix is the abstract base class for all classes representing the material properties of a
    dust medium. It implements the complete abstract interface defined by the MaterialMix base
    class for obtaining optical material properties. Depending on the scattering mode returned by
    the scatteringMode() function (which must be overridden by subclasses that support a scattering
    mode other than HenyeyGreenstein), the DustMix setup machinery requests the relevant optical
    dust properties from the subclass, and caches this information (possibly with some additional
    precalculated data) for later use. During the simulation, these properties can then be served
    up without further access to the subclass.

    In all cases, the properties served through the MaterialMix base class interface correspond to
    the properties of a single grain population that is representative for the complete dust mix.
    For dust mixes described by multiple grain populations, the optical properties should be
    integrated over the grain size distribution and accumulated across all grain populations. In
    the context of tracking photon paths through a dusty medium, using these integrated absorption
    and scattering cross sections and Mueller matrix coefficients is mathematically exact. In other
    words, for scattering modes MaterialPhaseFunction and SphericalPolarization, the representative
    grain approach does not involve an approximation. However, the calculation of a representative
    scattering asymmetry parameter \f$g\f$ for use with the Henyey-Greenstein scattering mode does
    involve a non-exact averaging procedure. Because the Henyey-Greenstein scattering phase
    function is non-physical to begin with, using a single approximate \f$g\f$ value for the
    complete dust mix is usually considered to be acceptable. For more information on the
    calculation of representative grain properties in SKIRT, refer to the description of the
    MultiGrainDustMix class.

    The representative grain properties offered by this class are insufficient to accurately
    calculate dust emission spectra for the dust mixture. This is so because the emission spectrum
    is a nonlinear function of (among many other things) the grain size, and thus a single grain
    cannot accurately represent a population with a (potentialy large) range of grain sizes. Again,
    refer to the description of the MultiGrainDustMix class for more information.

    The implementation of polarization by scattering in this class is based on the analysis
    presented by Peest at al. 2017 (A&A, 601, A92).

    <b>Extreme forward (or backward) Henyey-Greenstein scattering</b>

    In the X-ray wavelength range, dust grains exhibit extreme forward scattering with asymmetry
    parameter \f$g\f$ values very close to unity. Calculating the value of and sampling from the
    Henyey-Greenstein (HG) phase function becomes numerically unstable for \f$|g| > 1-10^{-6}\f$.
    Therefore, the implementation clips any larger \f$|g|\f$ values to that limit.

    More annoyingly, the bias weights for scattering peel-off photon packets can become very large,
    causing unacceptably high noise levels in the fluxes recorded by instruments. Consider the HG
    phase function, \f[ \Phi_\mathrm{HG}(\cos\theta) = \frac{1-g^2}{(1+g^2-2g\cos\theta)^{3/2}}.
    \f] Given the angle \f$\theta\f$ between the direction of the incoming photon packet and a
    given SKIRT instrument, the bias weight for a peel-off photon packet sent to the instrument is
    given by \f$w=\Phi_\mathrm{HG}(\cos\theta)\f$. For \f$g>0\f$, the largest bias weights are
    associated with peel-off photon packets emitted in the original propagation direction
    (\f$\theta=0\f$). Experiments show that the weigths should be kept under \f$\approx 10^3\f$ to
    achieve acceptable noise levels with a practical number of photon packets. Incidentally, for
    \f$g = 0.95\f$ we find \f$w\approx 780\f$, indicating the maximum asymmetry parameter value we
    could support given this criterion. For \f$g = 1-10^{-6}\f$, we find \f$w\approx 2 \times
    10^{12}\f$.

    To address this situation, we average the peel-off bias weight for a photon packet sent to an
    instrument over a portion of the phase function, as opposed to taking the value at a given
    angle. This introduces a certain amount of blur in recorded fluxes, while simultaneously
    decreasing the level of noise. The effect is reminiscent of that of an actual instrument's
    point spread function.

    We introduce the notation \f$\Psi(\alpha,\beta)\f$ for the definite integral of the HG phase
    function between two angles \f$0 \le \alpha \le \beta \le \pi\f$, and assuming \f$0 < |g| <
    1\f$, \f[ \Psi(\alpha,\beta) \equiv \int_{\cos\beta}^{\cos\alpha} \Phi_\mathrm{HG}(\cos\theta)
    \,\mathrm{d}\cos\theta = \frac{(1-g^2)}{g}\,\frac{(t_\beta-t_\alpha)}{t_\beta\,t_\alpha},\f]
    with \f[ t_x = \sqrt{1+g^2-2g\cos x}, \;\;x=\alpha,\beta. \f]

    We similarly use \f$\ell(\alpha,\beta)\f$ for the length of the interval between the cosines of
    two angles \f$0 \le \alpha \le \beta \le \pi\f$, \f[ \ell(\alpha,\beta) \equiv
    \cos\alpha-\cos\beta. \f]

    We further denote the half-opening angle of the instrument's solid angle as \f$\delta\f$, with
    \f$0 < \delta \ll \pi/2\f$. The averaged peel-off bias weight can then be written as, \f[
    \left<w\right> = \begin{cases} \displaystyle \frac{\Psi(\theta-\delta,\theta+\delta)}
    {\ell(\theta-\delta,\theta+\delta)} & \mathrm{for}\; \delta \le \theta \le \pi-\delta, \\[15pt]
    \displaystyle \frac{\Psi(0,\delta-\theta) + \Psi(0,\delta+\theta)} {\ell(0,\delta-\theta) +
    \ell(0,\delta+\theta)} & \mathrm{for}\; \theta < \delta, \\[15pt] \displaystyle
    \frac{\Psi(\theta-\delta,\pi) + \Psi(2\pi-\theta-\delta,\pi)} {\ell(\theta-\delta,\pi) +
    \ell(2\pi-\theta-\delta,\pi)} & \mathrm{for}\; \theta > \pi-\delta. \end{cases} \f]

    The last two expressions are needed to handle the cases where the averaging interval straddles
    zero or \f$\pi\f$ because the definitions above are not valid for angles outside that range.

    Investigation of these expressions shows that a half-opening angle of \f$\delta = 4^\circ\f$
    dampens the maximum averaged bias weight for \em all \f$g\f$ to \f$\left<w\right>\approx
    820\f$. For \f$g = 0.999\f$, the maximum weight reaches \f$\left<w\right>\approx 810\f$ and
    remains saturated from there on. In conclusion, the implementation uses the regular bias weight
    for \f$|g| \le 0.95\f$ and the averaged bias weight for larger values, with \f$\delta =
    4^\circ\f$. As a result, recorded fluxes will be blurred accordingly. */
class DustMix : public MaterialMix
{
    ITEM_ABSTRACT(DustMix, MaterialMix, "a dust mix")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function obtains and/or precalculates information used to serve optical properties for
        the dust mix to the simulation.

        First, it obtains the optical properties of the dust mix from the subclass on a wavelength
        grid with a range spanning all wavelengths in the simulation, and with sufficient
        resolution so that no further interpolation is needed. This includes the absorption and
        scattering cross sections, the scattering asymmetry parameter, and/or the Mueller matrix
        coefficients as required by the implemented scattering mode.

        If the simulation wavelength range extends beyond 10 cm, it is cut off at that value and
        any dust cross sections beyond 10 cm are forced to zero. None of the currently provided
        dust mixes offers optical properties beyond 10 cm, and historically any values outside the
        supported range are clamped to the nearest available value. This leads to substantially
        overestimated dust extinction in the radio wavelength range. Hence this "hack".

        Furthermore, if the simulation tracks the radiation field, this function precalculates the
        Planck-integrated absorption cross sections on an appropriate temperature grid. This
        information is used to obtain the equilibrium temperature of the material mix (or rather,
        of its representative grain population) in a given embedding radiation field. */
    void setupSelfAfter() override;

    /** This function must be implemented in each subclass to obtain the representative grain
        optical properties for the dust mix. Although the function arguments provide for all
        properties that may be required, the list of actually required properties for a particular
        subclass depends on the scattering mode advertised by the scatteringMode() function.

        The first two arguments respectively specify the wavelength grid and (if applicable) the
        scattering angle grid on which the properties must be tabulated. The function must store
        the absorption and scattering cross sections and the asymmetry parameter into the three
        subsequent output arrays, and the Mueller coefficients into the four subsequent output
        tables (indexed on wavelength and angle in that order). These arrays and tables will
        already have the appropriate size (corresponding to the input wavelength grids) when the
        function gets called. Finally, the function must return the dust mass per hydrogen atom for
        the dust mix.

        For the HenyeyGreenstein scattering mode, the function must output the cross sections and
        the asymmetry parameter, and it must leave the Mueller coeficient tables untouched. For
        scattering modes other than HenyeyGreenstein, the asymmetry parameter output array may be
        filled with "relevant" values (not used for calculations but possibly exposed to the user
        through a Probe), or it may be left untouched. For the SphericalPolarization scattering
        mode, the function must fill all four Mueller coeficient tables. For the
        MaterialPhaseFunction scattering mode, the function must fill only the first table and it
        must leave the other tables untouched.

        For the SpheroidalPolarization mode, the function must also fill the sigmaabsvv and
        sigmaabspolvv tables. */
    virtual double getOpticalProperties(const Array& lambdav, const Array& thetav, Array& sigmaabsv, Array& sigmascav,
                                        Array& asymmparv, Table<2>& S11vv, Table<2>& S12vv, Table<2>& S33vv,
                                        Table<2>& S34vv, ArrayTable<2>& sigmaabsvv, ArrayTable<2>& sigmaabspolvv) = 0;

    /** This function can be implemented in a subclass to initialize dust properties that are
        required to offer additional functionality. The argument specifies the wavelength grid on
        which the properties may be tabulated (i.e. the same grid as passed to the
        getOpticalProperties() function.

        The function is called by the DustMix class during setup after the getOpticalProperties()
        function has been called and its results have been processed. The function returns the
        number of bytes allocated by the subclass to support the extra features (this number is
        used for logging purposes). The default implementation of this function does nothing and
        returns zero. */
    virtual size_t initializeExtraProperties(const Array& lambdav);

    /** This function logs a warning message if the given range is smaller than the simulation
        wavelength range. It can (but does not have to) be called from subclasses that support dust
        properties for a limited wavelength range. */
    void informAvailableWavelengthRange(Range available);

    //======== Private support functions =======

private:
    /** This function returns the index in the private wavelength grid corresponding to the
        specified wavelength. The parameters for converting a wavelength to the appropriate index
        are stored in data members during setup. */
    int indexForLambda(double lambda) const;

    /** This function returns the index in the private scattering angle grid corresponding to the
        specified scattering angle. The parameters for converting a scattering angle to the
        appropriate index are built-in constants. */
    int indexForTheta(double theta) const;

    //======== Material type =======

public:
    /** This function returns the fundamental material type represented by this material mix, in
        other words it returns MaterialType::Dust. */
    MaterialType materialType() const override;

    /** This function returns true for this class because all dust mixes support secondary
        continuum emission. */
    bool hasContinuumEmission() const override;

    //======== Medium state setup =======

public:
    /** This function returns a list of StateVariable objects describing the specific state
        variables used by the receiving material mix. See the description of the
        MaterialMix::specificStateVariableInfo() function for more information.

        Dust mixes require just the standard specific state variable of type numberDensity , so
        this function returns a list containing a single item. */
    vector<StateVariable> specificStateVariableInfo() const override;

    //======== Low-level material properties =======

public:
    /** This function returns the dust mass \f$\mu\f$ per hydrogen atom for this dust mix. It
        returns a value that was pre-computed during setup. */
    double mass() const override;

    /** This function returns the absorption cross section per entity
        \f$\varsigma^{\text{abs}}_{\lambda}\f$ at wavelength \f$\lambda\f$. It retrieves the
        requested value from a table that was pre-computed during setup. */
    double sectionAbs(double lambda) const override;

    /** This function returns the scattering cross section per entity
        \f$\varsigma^{\text{sca}}_{\lambda}\f$ at wavelength \f$\lambda\f$. It retrieves the
        requested value from a table that was pre-computed during setup. */
    double sectionSca(double lambda) const override;

    /** This function returns the total extinction cross section per entity
        \f$\varsigma^{\text{ext}}_{\lambda} = \varsigma^{\text{abs}}_{\lambda} +
        \varsigma^{\text{sca}}_{\lambda}\f$ at wavelength \f$\lambda\f$. It retrieves the requested
        value from a table that was pre-computed during setup. */
    double sectionExt(double lambda) const override;

    /** This function returns the scattering asymmetry parameter \f$g_\lambda =
        \left<\cos\theta\right>\f$ at wavelength \f$\lambda\f$, which is used with the
        HenyeyGreenstein scattering mode. It retrieves the requested value from a table that was
        pre-computed during setup. */
    double asymmpar(double lambda) const override;

    //======== High-level photon life cycle =======

public:
    /** This function returns the absorption opacity \f$k^\text{abs}=n\varsigma^\text{abs}\f$ for
        the given wavelength and material state. The photon properties are not used. */
    double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the scattering opacity \f$k^\text{sca}=n\varsigma^\text{sca}\f$ for
        the given wavelength and material state. The photon properties are not used. */
    double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the extinction opacity \f$k^\text{ext}=k^\text{abs}+k^\text{sca}\f$
        for the given wavelength and material state. The photon properties are not used. */
    double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function calculates the contribution of the medium component associated with this
        material mix to the peel-off photon luminosity, polarization state, and wavelength shift
        for the given wavelength, geometry, material state, and photon properties. See the
        description of the MaterialMix::peeloffScattering() function for more information.

        For dust mixes, evaluation of the phase function depends on the scattering mode supported
        by supported by the dust mix, as defined by each subclass. For the most basic mode, the
        material mix provides a value for the scattering asymmetry parameter
        \f$g=\left<\cos\theta\right>\f$. A value of \f$g=0\f$ corresponds to isotropic scattering.
        Other values \f$-1\le g\le 1\f$ are substituted in the Henyey-Greenstein phase function,
        \f[ \Phi(\cos\theta) = \frac{1-g^2} {(1+g^2-2g\cos\theta)^{3/2}}. \f] For other scattering
        modes, the phase function provided by the material mix is invoked instead.

        In case polarization is supported in the current simulation configuration, the polarization
        state of the peel off photon packet is adjusted as well. The adjusted Stokes vector for a
        particular medium component is obtained as follows. The function rotates the Stokes vector
        from the reference direction in the previous scattering plane into the peel-off scattering
        plane, applies the Mueller matrix on the Stokes vector, and further rotates the Stokes
        vector from the reference direction in the peel-off scattering plane to the x-axis of the
        instrument to which the peel-off photon packet is headed. */
    void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs, Direction bfky,
                           const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function performs a scattering event on the specified photon packet in the spatial
        cell and medium component represented by the specified material state and the receiving
        material mix.

        For dust mixes, the operation depends on the scattering mode supported by the dust mix, as
        defined by each subclass. For the most basic mode, the dust mix provides a value for the
        scattering asymmetry parameter \f$g=\left<\cos\theta\right>\f$. For the value \f$g=0\f$,
        corresponding to isotropic scattering, a new direction is generated uniformly on the unit
        sphere. For other values \f$-1\le g\le 1\f$, a scattering angle \f$\theta\f$ is sampled
        from the Henyey-Greenstein phase function, \f[ \Phi(\cos\theta) =
        \frac{1-g^2}{(1+g^2-2g\cos\theta)^{3/2}}. \f] This can be accomplished as follows.
        Substituting \f$\mu=\cos\theta\f$, the probability distribution for \f$\mu\f$ (normalized
        to unity) becomes \f[ p(\mu)\,\text{d}\mu = \frac{1}{2} \,
        \frac{1-g^2}{(1+g^2-2g\mu)^{3/2}} \,\text{d}\mu \qquad -1\leq\mu\leq1 \f] We can use the
        transformation method to sample from this distribution. Given a uniform deviate
        \f$\mathcal{X}\f$, we need to solve \f[ {\mathcal{X}} = \int_{-1}^\mu p(\mu')\,\text{d}\mu'
        \f] Performing the integration and solving for \f$\mu\f$ yields \f[ \cos\theta = \mu =
        \frac{1+g^2-f^2}{2g} \quad\text{with}\quad f=\frac{1-g^2}{1-g+2g {\mathcal{X}}}
        \qquad\text{for}\; g\neq 0 \f] For other scattering modes, a function provided by the dust
        mix is invoked instead to obtain a random scattering direction for the photon packet.

        In case polarization is supported in the current simulation configuration, the polarization
        state of the photon packet is adjusted as well. Note that all media must either support
        polarization or not support it, mixing these support levels is not allowed. Compliance with
        this requirement is verified during setup of the simulation. The adjusted Stokes vector is
        obtained as follows, again using the randomly selected medium component. After obtaining
        the sampled scattering angles \f$\theta\f$ and \f$\phi\f$ from the material mix, the Stokes
        vector of the photon packet is rotated into the scattering plane and transformed by
        applying the Mueller matrix. Finally, the new direction is computed from the previously
        sampled \f$\theta\f$ and \f$\phi\f$ angles. */
    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======== Scattering implementation for dust mixes =======

public:
    /** This enumeration lists the scattering modes supported by the dust mix hierarchy. */
    enum class ScatteringMode {
        HenyeyGreenstein,
        MaterialPhaseFunction,
        SphericalPolarization,
        SpheroidalPolarization
    };

    /** This function returns the scattering mode supported by this dust mix. In the current
        implementation, this can be one of the following modes:

        - HenyeyGreenstein: the value returned by the asymmpar() function serves as the assymmetry
        parameter \f$g\f$ for the Henyey-Greenstein phase function. For a value of \f$g=0\f$, the
        Henyey-Greenstein phase function describes isotropic scattering.

        - MaterialPhaseFunction: this material type implements a custom phase function that depends
        only on the cosine of the scattering angle, for unpolarized radiation. Specifically, the
        phaseFunctionValueForCosine() and generateCosineFromPhaseFunction() functions are used to
        obtain the value of the phase function and to sample a scattering angle from it.

        - SphericalPolarization: this material type supports polarization through scattering by
        spherical particles. In this mode, the phase function depends on the polarization state of
        the incoming radiation, and the polarization state of the outgoing radiation must be
        updated appropriately. The phaseFunctionValue() and generateAnglesFromPhaseFunction()
        functions are used to obtain the value of the phase function and to sample a scattering
        angle from it, and the applyMueller() function is used to updated the polarization state.

        - SpheroidalPolarization: this material type supports polarization through scattering,
        absorption and emission by nonspherical, spheroidal particles. Currently, only \em emission
        is implemented and all other areas of the code treat spheroidal particles as if they were
        spherical.

        The implementation of this function in this base class returns the HenyeyGreenstein
        scattering mode as a default value. Subclasses that support another scattering mode must
        override this function and return the appropriate value. */
    virtual ScatteringMode scatteringMode() const;

private:
    /** This function is used with the MaterialPhaseFunction scattering mode, which assumes that
        the scattering phase function depends only on the cosine of the scattering angle. The
        function returns the value of the scattering phase function \f$\Phi_\lambda(\cos\theta)\f$
        at wavelength \f$\lambda\f$ for the specified scattering angle cosine \f$\cos\theta\f$,
        where the phase function is normalized as \f[\int_{-1}^1 \Phi_\lambda(\cos\theta)
        \,\mathrm{d}\cos\theta =2.\f]

        The phase function for unpolarized radiation is obtained from the general polarized case
        described for the phaseFunctionValue() function by setting the linear polarization degree
        \f$P_\text{L}\f$ to zero. The simplified function uses only the first Mueller matrix
        coefficient \f$S_{11}\f$. */
    double phaseFunctionValueForCosine(double lambda, double costheta) const;

    /** This function is used with the MaterialPhaseFunction scattering mode, which assumes that
        the scattering phase function depends only on the cosine of the scattering angle. The
        function generates a random scattering angle cosine sampled from the phase function
        \f$\Phi_\lambda(\cos\theta)\f$ at wavelength \f$\lambda\f$.

        The phase function for unpolarized radiation is obtained from the general polarized case
        described for the phaseFunctionValue() function by setting the linear polarization degree
        \f$P_\text{L}\f$ to zero. The simplified function no longer depends on \f$\phi\f$, and uses
        only the first Mueller matrix coefficient \f$S_{11}\f$. To sample the scattering angle
        \f$\theta\f$ from this phase function, we follow the same procedure as described for the
        generateAnglesFromPhaseFunction() function, with \f$P_\text{L}=0\f$. */
    double generateCosineFromPhaseFunction(double lambda) const;

    /** This function is used with the SphericalPolarization scattering mode. It returns the value
        of the scattering phase function \f$\Phi_\lambda(\theta,\phi)\f$ at wavelength
        \f$\lambda\f$ for the specified scattering angles \f$\theta\f$ and \f$\phi\f$, and for the
        specified incoming polarization state. The phase function is normalized as
        \f[\int\Phi_\lambda(\theta,\phi) \,\mathrm{d}\Omega =4\pi.\f]

        The phase function for scattering by spherical grains can be written as \f[
        \Phi(\theta,\phi) = N\,S_{11} \left( 1 + P_{\text{L}}\,\frac{S_{12}}{S_{11}}\cos2(\phi -
        \gamma) \right) \f] with \f[ N=\frac{1}{2\pi\int_0^\pi S_{11}\sin\theta\, \text{d}\theta},
        \f] where \f$P_\text{L}\f$ is the linear polarization degree and \f$\gamma\f$ the
        polarization angle of the incoming photon, and where the Mueller matrix coefficients
        \f$S_{xx}\f$ depend on both the photon wavelength \f$\lambda\f$ and the scattering angle
        \f$\theta\f$. */
    double phaseFunctionValue(double lambda, double theta, double phi, const StokesVector* sv) const;

    /** This function is used with the SphericalPolarization scattering mode. It generates random
        scattering angles \f$\theta\f$ and \f$\phi\f$ sampled from the phase function
        \f$\Phi_\lambda(\theta,\phi)\f$ at wavelength \f$\lambda\f$, and for the specified incoming
        polarization state. The results are returned as a pair of numbers in the order \f$\theta\f$
        and \f$\phi\f$.

        For scattering by spherical grains, we sample from the phase function listed for the
        phaseFunctionValue() function using the conditional probability technique. We reduce the
        phase function to the marginal distribution \f$\Phi(\theta)\f$, \f[ \Phi(\theta)
        =\int_0^{2\pi} \Phi(\theta,\phi)\,\text{d}\phi =2\pi\ N\,S_{11} =\frac{S_{11}}{\int_0^\pi
        S_{11}\sin\theta'\, \text{d}\theta'} \ . \f] We sample a random \f$\theta\f$ value from
        this distribution through numerical inversion, that is to say, by solving the equation
        \f[\label{eq:numInvTheta} {\cal{X}} =\frac{\int_0^\theta
        S_{11}\sin\theta'\,\text{d}\theta'}{\int_0^\pi S_{11}\sin\theta'\, \text{d}\theta'} \f] for
        \f$\theta\f$, where \f${\cal{X}}\f$ is a uniform deviate.

        Once we have selected a random scattering angle \f$\theta\f$, we sample a random azimuthal
        angle \f$\phi\f$ from the normalized conditional distribution, \f[ \Phi_\theta(\phi)
        =\frac{\Phi(\theta,\phi)}{\int_0^{2\pi} \Phi(\theta,\phi')\,\text{d}\phi'}
        =\frac{1}{2\pi}\left(1+ P_{\text{L}}\,\frac{S_{12}}{S_{11}}\cos 2(\phi - \gamma)\right),
        \f] where the ratio \f$S_{12}/S_{11}\f$ depends on the scattering angle \f$\theta\f$ in
        addition to wavelength.

        This can again be done through numerical inversion, by solving the equation \f[ {\cal{X}}
        =\int_{0}^{\phi}\Phi_{\theta}(\phi')\,\text{d}\phi' =\frac{1}{2\pi} \left( \phi +
        P_{\text{L}}\,\frac{S_{12}}{S_{11}} \sin\phi \cos(\phi - 2\gamma)\right) \f] for
        \f$\phi\f$, with \f${\cal{X}}\f$ being a new uniform deviate. */
    std::pair<double, double> generateAnglesFromPhaseFunction(double lambda, const StokesVector* sv) const;

    /** This function is used with the SphericalPolarization scattering mode. It applies the
        Mueller matrix transformation for the specified wavelength \f$\lambda\f$ and scattering
        angle \f$\theta\f$ to the given polarization state (which serves as both input and output
        for the function).

        For scattering by spherical grains, the Mueller matrix has only four independent
        coefficients, namely \f$S_{11}\f$, \f$S_{12}\f$, \f$S_{33}\f$, and \f$S_{34}\f$, which
        depend on both the photon wavelength \f$\lambda\f$ and the scattering angle \f$\theta\f$.
        These coefficients are obtained from the tables pre-computed during setup. */
    void applyMueller(double lambda, double theta, StokesVector* sv) const;

    //======== Secondary emission =======

public:
    /** This function returns the wavelength grid on which dust emission is discretized, i.e. the
        wavelength grid returned by the Configuration::dustEmissionWLG() function. */
    DisjointWavelengthGrid* emissionWavelengthGrid() const override;

    /** This function returns the emissivity spectrum \f$\varepsilon_{\ell'}\f$ (radiated power per
        unit of solid angle and per hydrogen atom) of the dust mix (or rather of the representative
        grain population corresponding to the dust mix) when it would be embedded in a given
        radiation field, assuming that the dust grains are in local thermal equilibrium. The input
        and output arrays are discretized on the wavelength grids returned by the
        Configuration::radiationFieldWLG() and Configuration::dustEmissionWLG() functions,
        repectively.

        The equilibrium emissivity of a representative grain population in an embedding radiation
        field \f$J_\lambda\f$ can be calculated as \f[ \varepsilon_\lambda =
        \varsigma_{\lambda,b}^{\text{abs}}\, B_\lambda(T_\text{eq}) \f] with \f$\mu\f$ the total
        dust mass of the dust mix, \f$\varsigma_{\lambda}^{\text{abs}}\f$ the absorption cross
        section of the representative grain, and \f$T_\text{eq}\f$ the equilibrium temperature of
        that grain, obtained from the energy balance equation as described for the
        DustMix::indicativeTemperature() function.

        The behavior of this function is undefined if the simulation does not track the radiation
        field, because in that case setup does not calculate the information on which this function
        relies. */
    Array emissivity(const Array& Jv) const override;

    /** This function returns the emission spectrum (radiated power per unit of solid angle) in the
        spatial cell and medium component represented by the specified material state and the
        receiving material mix when it would be embedded in the specified radiation field. For a
        dust mix, it returns the result of the emissivity() function multiplied by the hydrogen
        number density retrieved from the material state. */
    Array emissionSpectrum(const MaterialState* state, const Array& Jv) const override;

    /** This function returns an indicative temperature of the material mix when it would be
        embedded in a given radiation field. For dust mixes, it returns the equilibrium temperature
        \f$T_{\text{eq}}\f$ of the dust mix (or rather of the representative grain population
        corresponding to the dust mix) when it would be embedded in the radiation field specified
        by the mean intensities \f$(J_\lambda)_\ell\f$, which must be discretized on the
        simulation's radiation field wavelength grid as returned by the
        Configuration::radiationFieldWLG() function. If the specified \em Jv array is empty
        (because the simulation does not track the radiation field), this function returns zero.

        The equilibrium temperature is obtained from the energy balance equation, \f[ \int_0^\infty
        \varsigma^\text{abs}(\lambda) \,J_\lambda(\lambda) \,\text{d}\lambda = \int_0^\infty
        \varsigma^\text{abs}(\lambda) \,B_\lambda(T_\text{eq},\lambda) \,\text{d}\lambda, \f] where
        the left-hand side is integrated over the radiation field wavelength grid, and the
        right-hand side is precalculated during setup for a range of temperatures through
        integration over a built-in wavelength grid. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //======== Spheroidal grains =======

public:
    /** This function is intended for use with the SpheroidalPolarization mode. It returns the grid
        used for discretizing quantities that are a function of the scattering/emission angle
        \f$\theta\f$. The same grid is returned by all material mixes that have
        SpheroidalPolarization mode. */
    const Array& thetaGrid() const override;

    /** This function is intended for use with the SpheroidalPolarization mode. It returns the
        absorption cross sections per entity \f$\varsigma ^{\text{abs}} _{\lambda} (\theta)\f$ at
        wavelength \f$\lambda\f$ as a function of the emission angle \f$\theta\f$, discretized on
        the grid returned by the thetaGrid() function. */
    const Array& sectionsAbs(double lambda) const override;

    /** This function is intended for use with the SpheroidalPolarization mode. It returns the
        linear polarization absorption cross sections per entity \f$\varsigma ^{\text{abspol}}
        _{\lambda} (\theta)\f$ at wavelength \f$\lambda\f$ as a function of the emission angle
        \f$\theta\f$, discretized on the grid returned by the thetaGrid() function. */
    const Array& sectionsAbspol(double lambda) const override;

    //======================== Data Members ========================

private:
    // all data members are precalculated in setupSelfAfter()

    // wavelength grid (shifted to the left of the actually sampled points to approximate rounding)
    Array _lambdav;   // indexed on ell
    Range _required;  // the required wavelength range, i.e. the range of _lambdav before it was shifted

    // scattering angle grid
    Array _thetav;  // indexed on t

    // basic optical properties
    double _mu{0.};
    Array _sigmaabsv;  // indexed on ell
    Array _sigmascav;  // indexed on ell
    Array _sigmaextv;  // indexed on ell
    Array _asymmparv;  // indexed on ell

    // Mueller matrix coefficients
    Table<2> _S11vv;  // indexed on ell,t
    Table<2> _S12vv;  // indexed on ell,t
    Table<2> _S33vv;  // indexed on ell,t
    Table<2> _S34vv;  // indexed on ell,t

    // precalculated discretizations of (functions of) the scattering angles
    ArrayTable<2> _thetaXvv;  // indexed on ell and t
    Array _pfnormv;           // indexed on ell
    Array _phiv;              // indexed on f
    Array _phi1v;             // indexed on f
    Array _phisv;             // indexed on f
    Array _phicv;             // indexed on f

    // precalculated discretizations for spheroidal grains as a function of the emission angle
    ArrayTable<2> _sigmaabsvv;     // indexed on ell and t
    ArrayTable<2> _sigmaabspolvv;  // indexed on ell and t

    // equilibrium temperature and emission calculator
    EquilibriumDustEmissionCalculator _calc;
};

////////////////////////////////////////////////////////////////////

#endif
