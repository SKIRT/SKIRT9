/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTMIX_HPP
#define DUSTMIX_HPP

#include "ArrayTable.hpp"
#include "EquilibriumDustEmissionCalculator.hpp"
#include "MaterialMix.hpp"
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
    presented by Peest at al. 2017 (A&A, 601, A92). */
class DustMix : public MaterialMix
{
    ITEM_ABSTRACT(DustMix, MaterialMix, "a dust mix")
        ATTRIBUTE_TYPE_INSERT(DustMix, "Dust")
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

    //======== Private support functions =======

protected:
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
        other words it returns MaterialType::Dust. See the documentation of the MaterialMix class
        for more information. */
    MaterialType materialType() const override;

    //======== Basic material properties =======

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

    /** This function returns the scattering albedo \f$\varpi_\lambda =
        \varsigma_{\lambda}^{\text{sca}} / \varsigma_{\lambda}^{\text{ext}} =
        \kappa_{\lambda}^{\text{sca}} / \kappa_{\lambda}^{\text{ext}}\f$ at wavelength
        \f$\lambda\f$. It retrieves the requested value from a table that was pre-computed during
        setup. */
    double albedo(double lambda) const override;

    /** This function returns the scattering asymmetry parameter \f$g_\lambda =
        \left<\cos\theta\right>\f$ at wavelength \f$\lambda\f$, which is used with the
        HenyeyGreenstein scattering mode. It retrieves the requested value from a table that was
        pre-computed during setup. */
    double asymmpar(double lambda) const override;

    //======== Scattering with material phase function =======

public:
    /** This function returns the value of the scattering phase function
        \f$\Phi_\lambda(\cos\theta)\f$ at wavelength \f$\lambda\f$ for the specified scattering
        angle cosine \f$\cos\theta\f$, where the phase function is normalized as \f[\int_{-1}^1
        \Phi_\lambda(\cos\theta) \,\mathrm{d}\cos\theta =2.\f]

        The phase function for unpolarized radiation is obtained from the general polarized case
        described for the phaseFunctionValue() function by setting the linear polarization degree
        \f$P_\text{L}\f$ to zero. The simplified function uses only the first Mueller matrix
        coefficient \f$S_{11}\f$. */
    double phaseFunctionValueForCosine(double lambda, double costheta) const override;

    /** This function generates a random scattering angle cosine sampled from the phase function
        \f$\Phi_\lambda(\cos\theta)\f$ at wavelength \f$\lambda\f$.

        The phase function for unpolarized radiation is obtained from the general polarized case
        described for the phaseFunctionValue() function by setting the linear polarization degree
        \f$P_\text{L}\f$ to zero. The simplified function no longer depends on \f$\phi\f$, and uses
        only the first Mueller matrix coefficient \f$S_{11}\f$. To sample the scattering angle
        \f$\theta\f$ from this phase function, we follow the same procedure as described for the
        generateAnglesFromPhaseFunction() function, with \f$P_\text{L}=0\f$. */
    double generateCosineFromPhaseFunction(double lambda) const override;

    //======== Polarization through scattering by spherical particles =======

public:
    /** This function returns the value of the scattering phase function
        \f$\Phi_\lambda(\theta,\phi)\f$ at wavelength \f$\lambda\f$ for the specified scattering
        angles \f$\theta\f$ and \f$\phi\f$, and for the specified incoming polarization state. The
        phase function is normalized as \f[\int\Phi_\lambda(\theta,\phi) \,\mathrm{d}\Omega
        =4\pi.\f]

        The phase function for scattering by spherical grains can be written as \f[
        \Phi(\theta,\phi) = N\,S_{11} \left( 1 + P_{\text{L}}\,\frac{S_{12}}{S_{11}}\cos2(\phi -
        \gamma) \right) \f] with \f[ N=\frac{1}{2\pi\int_0^\pi S_{11}\sin\theta\, \text{d}\theta},
        \f] where \f$P_\text{L}\f$ is the linear polarization degree and \f$\gamma\f$ the
        polarization angle of the incoming photon, and where the Mueller matrix coefficients
        \f$S_{xx}\f$ depend on both the photon wavelength \f$\lambda\f$ and the scattering angle
        \f$\theta\f$. */
    double phaseFunctionValue(double lambda, double theta, double phi, const StokesVector* sv) const override;

    /** This function generates random scattering angles \f$\theta\f$ and \f$\phi\f$ sampled from
        the phase function \f$\Phi_\lambda(\theta,\phi)\f$ at wavelength \f$\lambda\f$, and for the
        specified incoming polarization state. The results are returned as a pair of numbers in the
        order \f$\theta\f$ and \f$\phi\f$.

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
    std::pair<double, double> generateAnglesFromPhaseFunction(double lambda, const StokesVector* sv) const override;

    /** This function applies the Mueller matrix transformation for the specified wavelength
        \f$\lambda\f$ and scattering angle \f$\theta\f$ to the given polarization state (which
        serves as both input and output for the function).

        For scattering by spherical grains, the Mueller matrix has only four independent
        coefficients, namely \f$S_{11}\f$, \f$S_{12}\f$, \f$S_{33}\f$, and \f$S_{34}\f$, which
        depend on both the photon wavelength \f$\lambda\f$ and the scattering angle \f$\theta\f$.
        These coefficients are obtained from the tables pre-computed during setup. */
    void applyMueller(double lambda, double theta, StokesVector* sv) const override;

    //======== Polarization through scattering, absorption and emission by spheroidal particles =======

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

    //======== Temperature and emission =======

public:
    /** This function returns the equilibrium temperature \f$T_{\text{eq}}\f$ of the dust mix (or
        rather of the representative grain population corresponding to the dust mix) when it would
        be embedded in the radiation field specified by the mean intensities
        \f$(J_\lambda)_\ell\f$, which must be discretized on the simulation's radiation field
        wavelength grid as returned by the Configuration::radiationFieldWLG() function.

        The equilibrium temperature is obtained from the energy balance equation, \f[ \int_0^\infty
        \varsigma^\text{abs}(\lambda) \,J_\lambda(\lambda) \,\text{d}\lambda = \int_0^\infty
        \varsigma^\text{abs}(\lambda) \,B_\lambda(T_\text{eq},\lambda) \,\text{d}\lambda, \f] where
        the left-hand side is integrated over the radiation field wavelength grid, and the
        right-hand side is precalculated during setup for a range of temperatures through
        integration over a built-in wavelength grid.

        The behavior of this function is undefined if the simulation does not track the radiation
        field, because in that case setup does not calculate the information on which this function
        relies. */
    double equilibriumTemperature(const Array& Jv) const override;

    /** This function returns the emissivity spectrum per hydrogen atom \f$\varepsilon_{\ell'}\f$
        of the dust mix (or rather of the representative grain population corresponding to the dust
        mix) when it would be embedded in a given radiation field, assuming that the dust grains
        are in local thermal equilibrium. The input and output arrays are discretized on the
        wavelength grids returned by the Configuration::radiationFieldWLG() and
        Configuration::dustEmissionWLG() functions, repectively.

        The equilibrium emissivity of a representative grain population in an embedding radiation
        field \f$J_\lambda\f$ can be calculated as \f[ \varepsilon_\lambda =
        \varsigma_{\lambda,b}^{\text{abs}}\, B_\lambda(T_\text{eq}) \f] with \f$\mu\f$ the total
        dust mass of the dust mix, \f$\varsigma_{\lambda}^{\text{abs}}\f$ the absorption cross
        section of the representative grain, and \f$T_\text{eq}\f$ the equilibrium temperature of
        that grain, obtained from the energy balance equation as described for the
        equilibriumTemperature() function.

        The behavior of this function is undefined if the simulation does not track the radiation
        field, because in that case setup does not calculate the information on which this function
        relies. */
    Array emissivity(const Array& Jv) const override;

    //======================== Data Members ========================

private:
    // all data members are precalculated in setupSelfAfter()

    // wavelength grid (shifted to the left of the actually sampled points to approximate rounding)
    Array _lambdav;  // indexed on ell

    // scattering angle grid
    Array _thetav;  // indexed on t

    // basic optical properties
    double _mu{0.};
    Array _sigmaabsv;  // indexed on ell
    Array _sigmascav;  // indexed on ell
    Array _sigmaextv;  // indexed on ell
    Array _albedov;    // indexed on ell
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
