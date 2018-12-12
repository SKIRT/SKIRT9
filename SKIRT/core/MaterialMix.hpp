/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MATERIALMIX_HPP
#define MATERIALMIX_HPP

#include "SimulationItem.hpp"
#include "Array.hpp"
class Random;
class StokesVector;
class WavelengthGrid;

////////////////////////////////////////////////////////////////////

/** MaterialMix is the abstract base class for all classes representing the concrete material
    properties of a specific transfer medium. The MaterialMix class hierarchy allows fundamentally
    different material types (e.g. dust, electrons, and hydrogen-dominated gas) to be implemented
    as part of a single framework. All MaterialMix subclasses, regardless of material type, must
    implement the public interface offered by this base class. This provides the material
    properties required for tracing photon packets through a material of this type, in other words,
    for processing absorption and scattering. The set of material properties that may be needed to
    calculate secondary emission spectra (e.g. thermal emission from dust grains) and/or medium
    state updates (e.g. hydrogen ionization fraction) differs between fundamental material types,
    and is thus offered by implementing one or more seperate interfaces (called "XxxMixInterface").

    Instances of MaterialMix subclasses are immutable after setup has been completed, so the same
    instance can be reused in multiple contexts.

    The MediumSystem object in each simulation maintains (directly or indirectly) three key data
    items to describe the contents of each spatial cell for each Medium component. The first is the
    cell volume. The second is a pointer to a MaterialMix, defining the properties of the material.
    The third is a number density, defining the amount of material present in the cell per unit of
    volume. The kind of physical entity being counted by the number density and the conversion
    from number density to mass density depend on the type of material, as indicated in the table
    below.

    Material type | Entity counted | Mass conversion
    --------------|----------------|-----------------------------
    Dust          | hydrogen atom  | dust mass per hydrogen atom
    Electrons     | electron       | electron mass
    Gas           | hydrogen atom  | gas mass per hydrogen atom

    The following table lists the key cell properties traced by the simulation, the key material
    properties offered by the interface in this base class, and some properties that can be
    derived from these.

    <TABLE>
    <TR><TD><B>Symbol</B></TD>  <TD><B>Units</B></TD>  <TD><B>Description</B></TD></TR>
    <TR><TD>\f$\Delta s\f$</TD>  <TD>\f$m\f$</TD>  <TD>Distance along a path</TD></TR>
    <TR><TD>\f$V\f$</TD>  <TD>\f$\text{m}^3\f$</TD>  <TD>Volume</TD></TR>
    <TR><TD>\f$n\f$</TD>  <TD>\f$\#\,\text{m}^{-3}\f$</TD>  <TD>Number density (of entities)</TD></TR>
    <TR><TD>\f$\mu\f$</TD>  <TD>\f$\text{kg}\,\#^{-1}\f$</TD>  <TD>Mass per entity</TD></TR>
    <TR><TD>\f$\varsigma\f$</TD>  <TD>\f$\text{m}^2\,\#^{-1}\f$</TD>  <TD>Cross section per entity</TD></TR>
    <TR><TD>\f$\mathcal{N}=n\Delta s\f$</TD> <TD>\f$\#\,\text{m}^{-2}\f$</TD>  <TD>Number column density</TD></TR>
    <TR><TD>\f$N=nV\f$</TD>  <TD>\f$\#\f$</TD>  <TD>Number (of entities)</TD></TR>
    <TR><TD>\f$\rho=n\mu\f$</TD>  <TD>\f$\text{kg}\,\text{m}^{-3}\f$</TD>  <TD>Mass density</TD></TR>
    <TR><TD>\f$\Sigma=n\mu\Delta s\f$</TD> <TD>\f$\text{kg}\,\text{m}^{-2}\f$</TD>  <TD>Mass column density</TD></TR>
    <TR><TD>\f$M=n\mu V\f$</TD>  <TD>\f$\text{kg}\f$</TD>  <TD>Mass</TD></TR>
    <TR><TD>\f$\kappa=\varsigma/\mu\f$</TD>  <TD>\f$\text{m}^2\,\text{kg}^{-1}\f$</TD>  <TD>Mass coefficient</TD></TR>
    <TR><TD>\f$k=n\varsigma\f$</TD>  <TD>\f$\text{m}^{-1}\f$</TD>  <TD>Opacity</TD></TR>
    <TR><TD>\f$\tau=n\varsigma\Delta s\f$</TD>  <TD>\f$1\f$</TD>  <TD>Optical depth</TD></TR>
    </TABLE>

    This class offers an interface to obtain the following basic material properties:
    - the mass per entity \f$\mu\f$
    - the absorption cross section per entity \f$\varsigma^{\text{abs}}_{\lambda}\f$
    - the scattering cross section per entity \f$\varsigma^{\text{sca}}_{\lambda}\f$
    - the total extinction cross section per entity \f$\varsigma^{\text{ext}}_{\lambda}
         = \varsigma^{\text{abs}}_{\lambda} + \varsigma^{\text{sca}}_{\lambda}\f$
    - the scattering albedo \f$\varpi_\lambda
         = \varsigma_{\lambda}^{\text{sca}} / \varsigma_{\lambda}^{\text{ext}}\f$
    - scattering phase function properties; depending on the supported scattering mode this may be:
      - the assymmetry parameter \f$g\f$ for the Henyey-Greenstein phase function
      - a custom phase function \f$\Phi_\lambda(\cos\theta)\f$ that depends only on the cosine of
        the scattering angle \f$\theta\f$
      - a custom phase function \f$\Phi_\lambda(\theta,\phi)\f$ that depends on both scattering angles
        \f$\theta\f$ and \f$\phi\f$, and on the polarization state of the incoming radiation
    - the equilibrium temperature \f$T_{\text{eq}}\f$ for a given embedding radiation field,
      using average material properties and assuming local thermal equilibrium conditions
*/
class MaterialMix : public SimulationItem
{
    ITEM_ABSTRACT(MaterialMix, SimulationItem, "a material mix")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function caches the simulation's random generator for use by subclasses. */
    void setupSelfBefore() override;

    //======== Material type =======

public:
    /** This enumeration lists the fundamental material types supported by the MaterialMix class
        hierarchy. */
    enum class MaterialType { Dust, Electrons, Gas };

    /** This function returns the fundamental material type represented by this material mix. See
        the documentation of the MaterialMix class for more information. */
    virtual MaterialType materialType() const = 0;

    /** This convenience function returns true if the fundamental material type represented by this
        material mix is Dust, and false otherwise. */
    bool isDust() const { return materialType() == MaterialType::Dust; }

    /** This convenience function returns true if the fundamental material type represented by this
        material mix is Electrons, and false otherwise. */
    bool isElectrons() const { return materialType() == MaterialType::Electrons; }

    /** This convenience function returns true if the fundamental material type represented by this
        material mix is Gas, and false otherwise. */
    bool isGas() const { return materialType() == MaterialType::Gas; }

    //======== Functionality levels =======

public:
    /** This enumeration lists the possible scattering modes offered by the public material mix
        interface. */
    enum class ScatteringMode { HenyeyGreenstein, MaterialPhaseFunction, SphericalPolarization };

    /** This function returns the scattering mode supported by this material mix. In the current
        implementation, this can be one of the following modes:

        - HenyeyGreenstein: the value returned by the asymmpar() function serves as the assymmetry
          parameter \f$g\f$ for the Henyey-Greenstein phase function. For a value of \f$g=0\f$,
          isotropic scattering is implemented directly (rather than subsituting zero into the
          Henyey-Greenstein phase function).

        - MaterialPhaseFunction: this material type implements a custom phase function that depends
          only on the cosine of the scattering angle, for unpolarized radiation. Specifically, the
          phaseFunctionValueForCosine() and generateCosineFromPhaseFunction() functions are used to
          obtain the value of the phase function and to sample a scattering angle from it.

        - SphericalPolarization: this material type supports polarization through scattering by
          spherical particles. In this mode, the phase function depends on the polarization state
          of the incoming radiation, and the polarization state of the outgoing radiation must be
          updated appropriately. The phaseFunctionValue() and generateAnglesFromPhaseFunction()
          functions are used to obtain the value of the phase function and to sample a scattering
          angle from it, and the applyMueller() function is used to updated the polarization state.
    */
    virtual ScatteringMode scatteringMode() const = 0;

    /** This convenience function returns true if this material mix uses and supports polarized
        radiation, and false otherwise. In the current implementation, the function returns true
        only if the scattering mode is SphericalPolarization. */
    bool hasPolarization() const { return scatteringMode() == ScatteringMode::SphericalPolarization; }

    //======== Basic material properties =======

public:
    /** This function returns the mass per entity \f$\mu\f$ for this material. The table below
        indicates the precise meaning of this number depending on the type of material being
        represented.

        Material type | Interpretation of mass() return value
        --------------|---------------------------------------
        Dust          | dust mass per hydrogen atom
        Electrons     | electron mass
        Gas           | gas mass per hydrogen atom
        */
    virtual double mass() const = 0;

    /** This function returns the absorption cross section per entity
        \f$\varsigma^{\text{abs}}_{\lambda}\f$ at wavelength \f$\lambda\f$. */
    virtual double sectionAbs(double lambda) const = 0;

    /** This function returns the scattering cross section per entity
        \f$\varsigma^{\text{sca}}_{\lambda}\f$ at wavelength \f$\lambda\f$. */
    virtual double sectionSca(double lambda) const = 0;

    /** This function returns the total extinction cross section per entity
        \f$\varsigma^{\text{ext}}_{\lambda} = \varsigma^{\text{abs}}_{\lambda} +
        \varsigma^{\text{sca}}_{\lambda}\f$ at wavelength \f$\lambda\f$. */
    virtual double sectionExt(double lambda) const = 0;

    /** This function returns the scattering albedo \f$\varpi_\lambda =
        \varsigma_{\lambda}^{\text{sca}} / \varsigma_{\lambda}^{\text{ext}} =
        \kappa_{\lambda}^{\text{sca}} / \kappa_{\lambda}^{\text{ext}}\f$ at
        wavelength \f$\lambda\f$. */
    virtual double albedo(double lambda) const = 0;

    /** This function is used only with the HenyeyGreenstein scattering mode. It returns the
        scattering asymmetry parameter \f$g_\lambda = \left<\cos\theta\right>\f$ at wavelength
        \f$\lambda\f$. This value serves as a parameter for the Henyey-Greenstein phase function.
        For a value of \f$g=0\f$, isotropic scattering is implemented directly. The default
        implementation in this base class returns 0. */
    virtual double asymmpar(double lambda) const;

    //======== Scattering with material phase function =======

public:
    /** This function is used with the MaterialPhaseFunction scattering mode, which assumes that
        the scattering phase function depends only on the cosine of the scattering angle. The
        function returns the value of the scattering phase function \f$\Phi_\lambda(\cos\theta)\f$
        at wavelength \f$\lambda\f$ for the specified scattering angle cosine \f$\cos\theta\f$,
        where the phase function is normalized as \f[\int_{-1}^1 \Phi_\lambda(\cos\theta)
        \,\mathrm{d}\cos\theta =2.\f] The default implementation in this base class returns one,
        corresponding to isotropic scattering. */
    virtual double phaseFunctionValueForCosine(double lambda, double costheta) const;

    /** This function is used with the MaterialPhaseFunction scattering mode, which assumes that
        the scattering phase function depends only on the cosine of the scattering angle. The
        function generates a random scattering angle cosine sampled from the phase function
        \f$\Phi_\lambda(\cos\theta)\f$ at wavelength \f$\lambda\f$. The default implementation in
        this base class returns a value sampled uniformly over the interval [-1,1], corresponding
        to isotropic scattering. */
    virtual double generateCosineFromPhaseFunction(double lambda) const;

    //======== Polarization through scattering by spherical particles =======

public:
    /** This function is used with the SphericalPolarization scattering mode. It returns the value
        of the scattering phase function \f$\Phi_\lambda(\theta,\phi)\f$ at wavelength
        \f$\lambda\f$ for the specified scattering angles \f$\theta\f$ and \f$\phi\f$, and for the
        specified incoming polarization state. The phase function is normalized as
        \f[\int\Phi_\lambda(\theta,\phi) \,\mathrm{d}\Omega =4\pi.\f] The default implementation in
        this base class returns one, corresponding to isotropic scattering. */
   virtual double phaseFunctionValue(double lambda, double theta, double phi, const StokesVector* sv) const;

    /** This function is used with the SphericalPolarization scattering mode. It generates random
        scattering angles \f$\theta\f$ and \f$\phi\f$ sampled from the phase function
        \f$\Phi_\lambda(\theta,\phi)\f$ at wavelength \f$\lambda\f$, and for the specified incoming
        polarization state. The results are returned as a pair of numbers in the order \f$\theta\f$
        and \f$\phi\f$. The default implementation in this base class arbitrarily returns two zero
        values. */
    virtual std::pair<double,double> generateAnglesFromPhaseFunction(double lambda, const StokesVector* sv) const;

    /** This function is used with the SphericalPolarization scattering mode. It applies the
        Mueller matrix transformation for the specified wavelength \f$\lambda\f$ and scattering
        angle \f$\theta\f$ to the given polarization state (which serves as both input and output
        for the function). The default implementation in this base class does nothing. */
    virtual void applyMueller(double lambda, double theta, StokesVector* sv) const;

    //======== Equilibrium Temperature =======

    /** This function returns the equilibrium temperature \f$T_{\text{eq}}\f$ (assuming LTE
        conditions) of the material mix when it would be embedded in the radiation field specified
        by the mean intensities \f$(J_\lambda)_\ell\f$, which must be discretized on the
        simulation's radiation field wavelength grid as returned by the
        Configuration::radiationFieldWLG() function. The default implementation in this base class
        returns zero. */
    virtual double equilibriumTemperature(const Array& Jv) const;

    //======================== Other Functions =======================

protected:
    /** This function returns the simulation's random generator as a service to subclasses. */
    Random* random() const { return _random; }

    //======================== Data Members ========================

private:
    // data member initialized in setupSelfBefore
    Random* _random{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
