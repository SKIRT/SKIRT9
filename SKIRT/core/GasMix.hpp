/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GASMIX_HPP
#define GASMIX_HPP

#include "MaterialMix.hpp"

////////////////////////////////////////////////////////////////////

/** This class makes it possible to normalize the gas by mass or number density, through the \c
    MaterialMix framework. The optical properties are dummy implementations, since the properties
    of the gas cannot simply be scaled with a factor. The optical properties are different in every
    cell, so they have to be handled by MediumSystem. If there is ever an update where the dust
    distribution can change depending on conditions in a cell, then we might be able to create a
    more unified approach for the gas and the dust. */
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
    /** This function returns the gas mass per hydrogen atom i.e. the mass of H. */
    double mass() const override;

    /** This function returns the absorption cross section per entity
        \f$\varsigma^{\text{abs}}_{\lambda}\f$ at wavelength \f$\lambda\f$. It is zero, as the gas
        properties depend on the cell, and need to be calculated using the Gas module. */
    double sectionAbs(double lambda) const override;

    /** This function returns the scattering cross section per entity
        \f$\varsigma^{\text{sca}}_{\lambda}\f$ at wavelength \f$\lambda\f$. It is zero, as the gas
        properties depend on the cell, and need to be calculated using the Gas module. */
    double sectionSca(double lambda) const override;

    /** This function returns the total extinction cross section per entity
        \f$\varsigma^{\text{ext}}_{\lambda} = \varsigma^{\text{abs}}_{\lambda} +
        \varsigma^{\text{sca}}_{\lambda}\f$ at wavelength \f$\lambda\f$. It is zero, as the gas
        properties depend on the cell, and need to be calculated using the Gas module. */
    double sectionExt(double lambda) const override;

    /** This function returns the scattering albedo \f$\varpi_\lambda =
        \varsigma_{\lambda}^{\text{sca}} / \varsigma_{\lambda}^{\text{ext}} =
        \kappa_{\lambda}^{\text{sca}} / \kappa_{\lambda}^{\text{ext}}\f$ at wavelength
        \f$\lambda\f$. It is zero, as the gas properties depend on the cell, and need to be
        calculated using the Gas module. */
    double albedo(double lambda) const override;

    //======== Temperature and emission =======

public:
    /** This function returns the equilibrium temperature \f$T_{\text{eq}}\f$ (assuming LTE
        conditions) of the gas when it would be embedded in the radiation field specified by the
        mean intensities \f$(J_\lambda)_\ell\f$, which must be discretized on the simulation's
        radiation field wavelength grid as returned by the Configuration::radiationFieldWLG()
        function. It is 1., as the gas properties depend on the cell, and need to be calculated
        using the Gas module. */
    double equilibriumTemperature(const Array& Jv) const override;

    /** This function returns the emissivity spectrum \f$\varepsilon_{\ell'}\f$ of the gas when it
        would be embedded in the radiation field specified by the mean intensities
        \f$(J_\lambda)_\ell\f$. The input radiation field must be discretized on the simulation's
        radiation field wavelength grid as returned by the Configuration::radiationFieldWLG()
        function. The output emissivity spectrum is discretized on a wavelength grid that depends
        on the material type. It is zero, as the gas properties depend on the cell, and need to be
        calculated using the Gas module. */
    Array emissivity(const Array& Jv) const override;
};

////////////////////////////////////////////////////////////////////

#endif
