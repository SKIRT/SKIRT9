/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LYANEUTRALHYDROGENMATERIALMIX_HPP
#define LYANEUTRALHYDROGENMATERIALMIX_HPP

#include "MaterialMix.hpp"

////////////////////////////////////////////////////////////////////

/** The LyaNeutralHydrogenMaterialMix class describes the material properties related to
    Lyman-alpha line transfer for a population of neutral hydrogen atoms, including support for
    polarization by scattering.

    The spatial distributions for both the mass density and the temperature of the neutral hydrogen
    gas must be defined by the input model and are considered to be constant during the simulation.
    In this context, this material mix offers a configuration property to specify a default gas
    temperature that is used by geometric media as a fixed temperature across the spatial domain.
    Imported media offer an \em importTemperature flag that allows defining a different gas
    temperature for each particle or cell. In this case, the default dust temperature configured
    for the material mix is ignored.

    This material mix also offers a configuration property to enable or disable support for
    polarization. It returns one of the \c Lya or \c LyaPolarization scattering modes depending on
    the configured value for this property. The dust-oriented functions for retrieving material
    properties are not used for the Lyman-alpha-specific aspects of the photon cycle, but some of
    them are used during setup, for example to normalize the mass of a geometric medium component
    based on optical depth. To this end, the scattering cross section is calculated using the
    configured default temperature and the absorption cross section is zero. */
class LyaNeutralHydrogenMaterialMix : public MaterialMix
{
    ITEM_CONCRETE(LyaNeutralHydrogenMaterialMix, MaterialMix, "neutral hydrogen for Lyman-alpha line transfer")
        ATTRIBUTE_TYPE_DISPLAYED_IF(LyaNeutralHydrogenMaterialMix, "Lya")

        PROPERTY_DOUBLE(defaultTemperature, "the default temperature of the neutral hydrogen gas")
        ATTRIBUTE_MIN_VALUE(defaultTemperature, "[3")  // gas temperature must be above local Universe T_CMB
        ATTRIBUTE_MAX_VALUE(defaultTemperature, "1e6]")
        ATTRIBUTE_DEFAULT_VALUE(defaultTemperature, "1e4")
        ATTRIBUTE_DISPLAYED_IF(defaultTemperature, "Level2")

        PROPERTY_BOOL(includePolarization, "include support for polarization")
        ATTRIBUTE_DEFAULT_VALUE(includePolarization, "false")
        ATTRIBUTE_DISPLAYED_IF(includePolarization, "Level2")

    ITEM_END()

    //======== Functionality levels =======

public:
    /** This function returns the fundamental material type represented by this material mix, which
        is MaterialType::Gas. See the documentation of the MaterialMix class for more information.
        */
    MaterialType materialType() const override;

    /** This function returns the scattering mode supported by this material mix, which is
        ScatteringMode::Lya or ScatteringMode::LyaPolarization depending on the value of the \em
        includePolarization flag. */
    ScatteringMode scatteringMode() const override;

    //======== Basic material properties =======

public:
    /** This function returns the mass of neutral hydrogen atom. */
    double mass() const override;

    /** This function returns the Lyman-alpha absorption cross section per hydrogen atom
        \f$\varsigma^\text{abs}_\alpha(\lambda, T_\text{def})\f$ which is trivially zero for all
        wavelengths and temperatures. */
    double sectionAbs(double lambda) const override;

    /** This function returns the Lyman-alpha scattering cross section per hydrogen atom
        \f$\varsigma^\text{sca}_\alpha(\lambda, T_\text{def})\f$ at the given wavelength and using
        the default gas temperature configured for this material mix. */
    double sectionSca(double lambda) const override;

    /** This function returns the total Lyman-alpha extinction cross section per hydrogen atom
        \f$\varsigma^\text{ext}_\alpha(\lambda, T_\text{def})\f$ at the given wavelength and using
        the default gas temperature configured for this material mix. The extinction cross section
        is identical to the scattering cross section because the absorption cross section is zero.
        */
    double sectionExt(double lambda) const override;

    /** This function returns the Lyman-alpha scattering albedo \f$\varpi_\alpha(\lambda,
        T_\text{def})\f$ for the hydrogen population, which is trivially equal to one for all
        wavelengths and temperatures. */
    double albedo(double lambda) const override;

    //======== Temperature and emission =======

    /** This function returns the equilibrium temperature \f$T_{\text{eq}}\f$ of the material mix
        when it would be embedded in a given radiation field. Because the hydrogen atoms absorb nor
        emit at the Lyman-alpha line (at least in our treatment), the implementation of this
        function in this class ignores the input radiation field and always returns the default
        temperature configured for this material mix. This is in fact a hack to allow retrieval of
        this temperature without casting the material mix object. */
    double equilibriumTemperature(const Array& Jv) const override;

    /** This function returns the emissivity spectrum \f$\varepsilon_{\ell'}\f$ of the material mix
        when it would be embedded in a given radiation field. Because the hydrogen atoms absorb nor
        emit at the Lyman-alpha line (at least in our treatment), the implementation of this
        function in this class ignores the input radiation field and always returns an empty array.
        */
    Array emissivity(const Array& Jv) const override;
};

////////////////////////////////////////////////////////////////////

#endif
