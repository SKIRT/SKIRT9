/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEANTABULATEDDUSTMIX_HPP
#define MEANTABULATEDDUSTMIX_HPP

#include "DustMix.hpp"
#include "Array.hpp"

////////////////////////////////////////////////////////////////////

/** MeanTabulatedDustMix is an abstract class for representing basic dust mixes described by
    tabulated properties for a single representative grain using the Henyey-Greenstein scattering
    mode.

    Specifically, the tabulated properties include the extinction mass coefficient
    \f$\kappa^\text{ext}_\lambda\f$, the scattering albedo \f$\varpi_\lambda\f$ and the scattering
    asymmetry parameter \f$g_\lambda\f$, as a function of wavelength \f$\lambda\f$. Because a basic
    dust mix such as this one is usually used in isolation and the dust distribution is normalized
    to a given optical depth or total dust mass, the value of the extinction coefficient is
    essentially scale free. Still, the dust mass per hydrogen atom \f$\mu\f$ may be specified to
    set the absolute scale of the property, so that the cross sections listed by some of the probes
    have an appropriately scaled value.

    Property values outside of the tabulated wavelength range are clamped to the nearest border
    value. As a special-case consequence, if only a single wavelength is tabulated, the properties
    are considered to be constant for all wavelengths.

    The subclass must load the tabulated data, and this abstract class handles everything else. */
class MeanTabulatedDustMix : public DustMix
{
    ITEM_ABSTRACT(MeanTabulatedDustMix, MaterialMix, "a basic dust mix with properties tabulated by the user")
        ATTRIBUTE_TYPE_DISPLAYED_IF(MeanTabulatedDustMix, "Level2")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function asks the subclass to load the tabulated data. */
    void setupSelfBefore() override;

    /** This function must be implemented in each subclass to store the wavelengths and the
        corresponding tabulated properties in the array arguments, and to return the dust mass per
        hydrogen atom. The function must guarantee that all arrays have the same size. */
    virtual double getDustProperties(Array& lambdav, Array& kappaextv, Array& albedov, Array& asymmparv) const = 0;

    //======================== Other Functions =======================

public:
    /** This function returns the scattering mode supported by this material mix. For this class,
        it returns the HenyeyGreenstein mode. */
    ScatteringMode scatteringMode() const override;

    /** This function returns the dust mass per hydrogen atom \f$\mu\f$ for this dust mix. */
    double mass() const override;

    /** This function returns the absorption cross section per hydrogen atom
        \f$\varsigma^{\text{abs}}_{\lambda}\f$ of the dust mix at wavelength \f$\lambda\f$. */
    double sectionAbsSelf(double lambda) const override;

    /** This function returns the scattering cross section per hydrogen atom
        \f$\varsigma^{\text{sca}}_{\lambda}\f$ of the dust mix at wavelength \f$\lambda\f$. */
    double sectionScaSelf(double lambda) const override;

    /** This function returns the scattering asymmetry parameter \f$g_\lambda =
        \left<\cos\theta\right>\f$ at wavelength \f$\lambda\f$, or if this value is unkown, it
        returns zero (corresponding to isotropic scattering). */
    double asymmpar(double lambda) const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    size_t _num{0};
    Array _lambdav;
    Array _sectionAbsv;
    Array _sectionScav;
    Array _asymmparv;
    double _mu{0.};
};

////////////////////////////////////////////////////////////////////

#endif
