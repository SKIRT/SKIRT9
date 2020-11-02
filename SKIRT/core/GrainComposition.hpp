/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GRAINCOMPOSITION_HPP
#define GRAINCOMPOSITION_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** GrainComposition is an abstract class that represents the optical and calorimetric properties
    of a population of dust grains with a given material composition.

    The basic optical properties are provided for arbitrary grain sizes and at arbitrary
    wavelengths; in practice, they are defined on a two-dimensional grid of wavelengths
    \f$\lambda_k\f$ and grain sizes \f$a_i\f$. The basic optical properties include the absorption
    and scattering efficiencies \f$Q_{k,i}^{\text{abs}}\f$ and \f$Q_{k,i}^{\text{sca}}\f$ and the
    scattering phase function asymmetry parameter \f$g_{k,i}\f$. All of these quantities are
    dimensionless.

    Optionally, a grain composition may provide the four coefficients \f$S_{xx}\f$ of the Mueller
    matrix for spherical grains, as a function of wavelength, grain size, and scattering angle. A
    dust mix containing grain compositions that provide the Mueller matrix coefficients can support
    the MaterialPhaseFunction and SphericalPolarization scattering modes, in addition to the
    HenyeyGreenstein scattering mode (see the MaterialMix class for more information on scattering
    modes).

    Furthermore, a grain composition may optionally provide additional properties relevant for
    calculating the polarization effects of scattering, absorption and emission by spheroidal
    grains. The current implementation supports only the properties relevant for polarized
    emission, i.e the absorption efficiencies \f$Q^\mathrm{abs}(a,\lambda,\theta)\f$ and the
    corresponding linear polarization efficiencies \f$Q^\mathrm{abspol}(a,\lambda,\theta)\f$ as a
    function of grain size \f$a\f$, wavelength \f$\lambda\f$ and emission angle \f$\theta\f$. For
    consistency, the \f$Q^\mathrm{abs}\f$ values given here, when averaged over the emission angle,
    should match the basic \f$Q^\mathrm{abs}\f$ mentioned above.

    The calorimetric properties consist of the specific enthalpy (internal energy) per unit volume,
    given in J/m3, at arbitrary temperature. The specific enthalpy is obtained by integrating the
    specific heat capacity of the material over the temperature range, using an arbitrary zero
    point. In practice, the specific enthalpy is defined on a grid of temperatures \f$T_t\f$
    resulting in a set of values \f$h_t\f$.

    A final, key property is the bulk mass density \f$\rho_\text{bulk}\f$ of the grain material, a
    single value given in kg/m3. To obtain the specific enthalpy per unit of mass (J/kg) from the
    provided specific enthalpy per unit volume (J/m3), one needs to divide the provided values by
    the bulk mass density (in kg/m3). Also, the same bulk mass density value must be used for
    calculating the dust mass from the size distribution to obtain consistent results.

    The GrainComposition class provides a public interface for retrieving the resource names for
    the various tables mentioned above, and for retrieving the value of the bulk mass density.
    Subclasses are required to implement these functions appropriately. */
class GrainComposition : public SimulationItem
{
    ITEM_ABSTRACT(GrainComposition, SimulationItem, "a dust grain composition")
    ITEM_END()

public:
    /** This function returns a brief human-readable identifier for the type of grain composition
        represented by the instance. The identifier is \em not allowed to contain white space. */
    virtual string name() const = 0;

    /** This function returns the bulk mass density \f$\rho_\text{bulk}\f$ of the grain material.
        */
    virtual double bulkDensity() const = 0;

    /** This function returns the name of the stored table resource tabulating the basic optical
        properties (absorption and scattering efficiencies and asymmetry parameter) as a function
        of wavelength and grain size. The asymmetry parameter values are used only with the
        HenyeyGreenstein scattering mode. */
    virtual string resourceNameForOpticalProps() const = 0;

    /** This function returns the name of the stored table resource tabulating the coefficients of
        the Mueller matrix as a function of wavelength, grain size, and scattering angle. This
        function is invoked for the MaterialPhaseFunction scattering mode (to obtain \f$S_{11}\f$ )
        and for the SphericalPolarization scattering mode (to obtain all four \f$S_{xx}\f$
        coefficients). The default implementation in this base class returns the empty string,
        indicating that only the HenyeyGreenstein scattering mode can be used. */
    virtual string resourceNameForMuellerMatrix() const;

    /** This function returns information on the resources required for implementing thermal
        emission from aligned spheriodal grains. It is invoked for the SpheroidalPolarization
        scattering mode. If the grain composition does not support this mode, the function returns
        false and the output arguments remain unchanghed. If the grain composition does support
        this mode, the function returns true and the output arguments are updated as follows:

        - \em resource: true if the specified tables are built-in resources or false if they are
        provided as user input files.

        - \em interpol: the interpolation fraction between the two tables in range [0,1]; the value
        is used to linearly interpolate between the two tables. A value of 0 means that the second
        table is ignored, a value of 1 means that the first table is ignored. A value in between is
        used to linearly interpolate between the two tables.

        - \em tableName1: the name of the first stored table (resource or user input file); this
        name must always be provided (i.e. non-empty and valid) even if \em interpol is equal to 1.

        - \em tableName2: the name of the second stored table (resource or user input file); this
        name may be missing (i.e. the empty string) if \em interpol is equal to 0.

        Each of the stored tables tabulates the absorption efficiencies
        \f$Q^\mathrm{abs}(a,\lambda,\theta)\f$ and the corresponding linear polarization
        efficiencies \f$Q^\mathrm{abspol}(a,\lambda,\theta)\f$ as a function of grain size \f$a\f$,
        wavelength \f$\lambda\f$ and emission angle \f$\theta\f$.

        The default implementation in this base class returns false and does not change the output
        arguments, indicating that there is no support for spheroidal grains. */
    virtual bool resourcesForSpheroidalEmission(bool& resource, double& interpol, string& tableName1,
                                                string& tableName2) const;

    /** This function returns the name of the stored table resource tabulating the specific
        enthalpies per unit volume as a function of temperature. */
    virtual string resourceNameForEnthalpies() const = 0;
};

////////////////////////////////////////////////////////////////////

#endif
