/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SINGLEGRAINDUSTMIX_HPP
#define SINGLEGRAINDUSTMIX_HPP

#include "DustMix.hpp"

////////////////////////////////////////////////////////////////////

/** SingleGrainDustMix is an abstract class implementing a dust mix described by a single
    representative grain, with or without support for polarization by scattering. This base class
    includes the implementations of the required functions to retrieve the optical properties from
    stored table resources. Subclasses must merely provide the names of the relevant resource files
    and implement the scatteringMode() function if they support a scattering mode other than
    Henyey-Greenstein. */
class SingleGrainDustMix : public DustMix
{
    ITEM_ABSTRACT(SingleGrainDustMix, DustMix, "a dust mix described by a single representative grain")
        ATTRIBUTE_TYPE_ALLOWED_IF(SingleGrainDustMix, "!StochasticDustEmission")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function is invoked by the DustMix base class to obtain the representative grain
        optical properties for this dust mix. The first two arguments respectively specify the
        wavelength grid and (if applicable) the scattering angle grid on which the properties must
        be tabulated. The output arrays and tables will already have the appropriate size
        (corresponding to the input wavelength grids) when the function gets called.

        For this class, this function first invokes the resourceNameForOpticalProps() function
        provided by each subclass to obtain the appropriate resource name for the basic optical
        properties. It opens the associated stored table resource, retrieves the optical properties
        on the requested wavelength grid, and stores them into the corresponding output arrays.

        Subsequently, for scattering modes other than HenyeyGreenstein, the function invokes the
        resourceNameForMuellerMatrix() function provided by each subclass to obtain the appropriate
        resource name for the Mueller matrix coefficients. For the SphericalPolarization scattering
        mode, the function fills all four tables. For the MaterialPhaseFunction scattering mode,
        the function fills only the first table and leaves the other tables untouched.

        Finally, the function returns the dust mass per hydrogen atom for the dust mix. */
    double getOpticalProperties(const Array& lambdav, const Array& thetav, Array& sigmaabsv, Array& sigmascav,
                                Array& asymmparv, Table<2>& S11vv, Table<2>& S12vv, Table<2>& S33vv, Table<2>& S34vv,
                                ArrayTable<2>& sigmaabsvv, ArrayTable<2>& sigmaabspolvv) override;

    /** This function must be implemented in a subclass to return the name of the stored table
        resource tabulating the basic optical properties (cross sections and asymmetry parameter)
        as a function of wavelength. The asymmetry parameter values are used only with the
        HenyeyGreenstein scattering mode (although they might be exposed to the user by a Probe). */
    virtual string resourceNameForOpticalProps() const = 0;

    /** This function must be implemented in a subclass to return the name of the stored table
        resource tabulating the coefficients of the Mueller matrix as a function of wavelength and
        scattering angle. This function is invoked for the MaterialPhaseFunction scattering mode
        (to obtain \f$S_{11}\f$ ) and for the SphericalPolarization (to obtain all four
        \f$S_{xx}\f$ coefficients). The default implementation in this base class returns the empty
        string, which is acceptable only with the HenyeyGreenstein scattering mode. */
    virtual string resourceNameForMuellerMatrix() const;
};

////////////////////////////////////////////////////////////////////

#endif
