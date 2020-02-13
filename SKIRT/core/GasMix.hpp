/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GASMIX_HPP
#define GASMIX_HPP

////////////////////////////////////////////////////////////////////

/** This class ties the gas to the MaterialMix framework, so that its optical properties can be
    used in the photon loops. */
class GasMix : public MaterialMix
{
    ITEM_CONCRETE(GasMix, MaterialMix, "a gas mix - only one is allowed")
    ITEM_END()

    //======== Material type =======

public:
    /** This function returns the fundamental material type represented by this material mix, in
        other words it returns MaterialType::Gas. See the documentation of the MaterialMix class
        for more information. */
    MaterialType materialType() const override;

    //======== Basic material properties =======

public:
    /** This function returns the gas mass per hydrogen atom i.e. 1. */
    double mass() const override;

    /** This function returns the absorption cross section per entity
        \f$\varsigma^{\text{abs}}_{\lambda}\f$ at wavelength \f$\lambda\f$. */
    double sectionAbs(double lambda) const override;

    /** This function returns the scattering cross section per entity
        \f$\varsigma^{\text{sca}}_{\lambda}\f$ at wavelength \f$\lambda\f$. */
    double sectionSca(double lambda) const override;

    /** This function returns the total extinction cross section per entity
        \f$\varsigma^{\text{ext}}_{\lambda} = \varsigma^{\text{abs}}_{\lambda} +
        \varsigma^{\text{sca}}_{\lambda}\f$ at wavelength \f$\lambda\f$. */
    double sectionExt(double lambda) const override;

    /** This function returns the scattering albedo \f$\varpi_\lambda =
        \varsigma_{\lambda}^{\text{sca}} / \varsigma_{\lambda}^{\text{ext}} =
        \kappa_{\lambda}^{\text{sca}} / \kappa_{\lambda}^{\text{ext}}\f$ at
        wavelength \f$\lambda\f$. */
    double albedo(double lambda) const override;

    //======== Temperature and emission =======

public:
    /** This function returns the equilibrium temperature \f$T_{\text{eq}}\f$ (assuming LTE
        conditions) of the gas when it would be embedded in the radiation field specified by the
        mean intensities \f$(J_\lambda)_\ell\f$, which must be discretized on the simulation's
        radiation field wavelength grid as returned by the Configuration::radiationFieldWLG()
        function.

        The argument is ignored, as the temperature is calculated during a gas state update. */
    double equilibriumTemperature(const Array& Jv) const override;

    /** This function returns the emissivity spectrum \f$\varepsilon_{\ell'}\f$ of the gas
        when it would be embedded in the radiation field specified by the mean intensities
        \f$(J_\lambda)_\ell\f$. The input radiation field must be discretized on the simulation's
        radiation field wavelength grid as returned by the Configuration::radiationFieldWLG()
        function. The output emissivity spectrum is discretized on a wavelength grid that depends
        on the material type. For more information, refer to the documentation of this function for
        each material type.

        The argument is ignored (but might be used later), as the emissivity can be calculated from
        the current gas state. */
    Array emissivity(const Array& Jv) const override;
};

////////////////////////////////////////////////////////////////////

#endif
