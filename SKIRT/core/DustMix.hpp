/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTMIX_HPP
#define DUSTMIX_HPP

#include "MaterialMix.hpp"

////////////////////////////////////////////////////////////////////

/** DustMix is the abstract base class for all classes representing the material properties of a
    dust medium. Other than implementing the materialType() function, its main purpose is to
    precalculate frequently used optical dust properties on the appropriate grids, so that this
    information can be retrieved much faster by the main simulation code. */
class DustMix : public MaterialMix
{
    ITEM_ABSTRACT(DustMix, MaterialMix, "a dust mix")
        ATTRIBUTE_TYPE_INSERT(DustMix, "Dust")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function precalculates information used to accelerate certain operations. First, it
        calculates the optical properties of the dust mix on a wavelength grid with a range
        spanning all wavelengths in the simulation, and with sufficient resolution so that no
        further interpolation is needed. Furthermore, if the simulation tracks the radiation field,
        it precalculates the Planck-integrated absorption cross sections on an appropriate
        temperature grid. This information is used to obtain the equilibrium temperature of the
        material mix in a given embedding radiation field.

        This function invokes the "Self" versions of the cross section retrieval functions, which
        should be implemented in each subclass. */
    void setupSelfAfter() override;

    //======== Material type =======

public:
    /** This function returns the fundamental material type represented by this material mix, in
        other words it returns MaterialType::Dust. See the documentation of the MaterialMix class
        for more information. */
    MaterialType materialType() const override;

    //======== Basic material properties =======

public:
    /** This function returns the absorption cross section per entity
        \f$\varsigma^{\text{abs}}_{\lambda}\f$ at wavelength \f$\lambda\f$. It retrieves the
        requested value from a data array that was pre-computed during setup. */
    double sectionAbs(double lambda) const override;

    /** This function returns the scattering cross section per entity
        \f$\varsigma^{\text{sca}}_{\lambda}\f$ at wavelength \f$\lambda\f$. It retrieves the
        requested value from a data array that was pre-computed during setup. */
    double sectionSca(double lambda) const override;

    /** This function returns the total extinction cross section per entity
        \f$\varsigma^{\text{ext}}_{\lambda} = \varsigma^{\text{abs}}_{\lambda} +
        \varsigma^{\text{sca}}_{\lambda}\f$ at wavelength \f$\lambda\f$. It retrieves the requested
        value from a data array that was pre-computed during setup. */
    double sectionExt(double lambda) const override;

    /** This function returns the scattering albedo \f$\varpi_\lambda =
        \varsigma_{\lambda}^{\text{sca}} / \varsigma_{\lambda}^{\text{ext}} =
        \kappa_{\lambda}^{\text{sca}} / \kappa_{\lambda}^{\text{ext}}\f$ at wavelength
        \f$\lambda\f$. It retrieves the requested value from a data array that was pre-computed
        during setup. */
    double albedo(double lambda) const override;

protected:
    /** This function determines and returns the absorption cross section per entity
        \f$\varsigma^{\text{abs}}_{\lambda}\f$ at wavelength \f$\lambda\f$. It must be implemented
        in each subclass and is invoked during the setup procedure of this abstract class to
        precompute the cross sections on an appropriate wavelength grid. */
    virtual double sectionAbsSelf(double lambda) const = 0;

    /** This function determines and returns the scattering cross section per entity
        \f$\varsigma^{\text{sca}}_{\lambda}\f$ at wavelength \f$\lambda\f$. It must be implemented
        in each subclass and is invoked during the setup procedure of this abstract class to
        precompute the cross sections on an appropriate wavelength grid. */
    virtual double sectionScaSelf(double lambda) const = 0;

    //======== Equilibrium Temperature =======

public:
    /** This function returns the equilibrium temperature \f$T_{\text{eq}}\f$ of the material mix
        when it would be embedded in the radiation field specified by the mean intensities
        \f$(J_\lambda)_\ell\f$, which must be discretized on the simulation's radiation field
        wavelength grid as returned by the Configuration::radiationFieldWLG() function.

        The equilibrium temperature is obtained from the energy balance equation, \f[ \int_0^\infty
        \varsigma^\text{abs}(\lambda) \,J_\lambda(\lambda) \,\text{d}\lambda = \int_0^\infty
        \varsigma^\text{abs}(\lambda) \,B_\lambda(T_\text{eq},\lambda) \,\text{d}\lambda, \f] where
        the left-hand side is integrated over the radiation field wavelength grid, and the
        right-hand side is precalculated during setup for a range of temperatures through
        integration over a built-in wavelength grid.

        The behavior of this function is undefined if the simulation does not track the radiation
        field, because in that case setup does not calcalate the info on which this function
        relies. */
    double equilibriumTemperature(const Array& Jv) const override;

    //======================== Data Members ========================

private:
    // cross sections - precalculated in setupSelfAfter()
    double _logLambdaOffset{0.};
    double _logLambdaFactor{0.};
    int _maxLogLambda{0};
    Array _sectionAbs;
    Array _sectionSca;
    Array _sectionExt;
    Array _albedo;

    // equilibrium temperature info - precalculated in setupSelfAfter()
    WavelengthGrid* _radiationFieldWLG{nullptr};  // the radiation field wavelength grid
    Array _sigmaabsv;   // absorption cross sections on the above wavelength grid
    Array _Tv;          // temperature grid for the Planck-integrated absorption cross sections
    Array _planckabsv;  // Planck-integrated absorption cross sections for each of the temperatures in the above grid
};

////////////////////////////////////////////////////////////////////

#endif
