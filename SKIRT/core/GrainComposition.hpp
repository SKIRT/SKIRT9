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

    /** This function returns the name of the stored table resource tabulating the specific
        enthalpies per unit volume as a function of temperature. */
    virtual string resourceNameForEnthalpies() const = 0;
};

////////////////////////////////////////////////////////////////////

#endif
