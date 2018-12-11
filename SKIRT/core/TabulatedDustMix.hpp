/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TABULATEDDUSTMIX_HPP
#define TABULATEDDUSTMIX_HPP

#include "DustMix.hpp"
#include "Array.hpp"

////////////////////////////////////////////////////////////////////

/** TabulatedDustMix is an abstract class for representing basic dust mixes described by tabulated
    properties for a single representative grain using the Henyey-Greenstein scattering mode.

    Specifically, the tabulated properties include the extinction mass coefficient
    \f$\kappa^\text{ext}_\lambda\f$, the scattering albedo \f$\varpi_\lambda\f$ and the scattering
    asymmetry parameter \f$g_\lambda\f$, as a function of wavelength \f$\lambda\f$. Because a basic
    dust mix such as this one is usually used in isolation and the dust distribution is normalized
    to a given optical depth or total dust mass, the value of the extinction coefficient is
    essentially scale free. Still, the dust mass per hydrogen atom \f$\mu\f$ may be specified to
    set the absolute scale of the property, so that the cross sections listed by some of the probes
    have an appropriately scaled value.

    Property values outside of the tabulated wavelength range are considered to be zero.

    The subclass must load the tabulated data, and this abstract class handles everything else. */
class TabulatedDustMix : public DustMix
{
    ITEM_ABSTRACT(TabulatedDustMix, MaterialMix, "a basic dust mix with properties tabulated by the user")
        ATTRIBUTE_TYPE_DISPLAYED_IF(TabulatedDustMix, "Level2")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function is invoked by the DustMix base class to obtain the representative grain
        optical properties for this dust mix. The first argument specifies the wavelength grid on
        which the absorption and scattering cross sections and the asymmetry parameter must be
        tabulated. The second argument is not used. The function store these tabulated properties
        into the three subsequent output arrays, which will already have the same length as the
        input wavelength grid when the function gets called. The Mueller matrix tables remain
        untouched. Finally, the function returns the dust mass per hydrogen atom for the dust mix.

        This function in turn invokes the getDustProperties() function that must be implemented by
        each TabulatedDustMix subclass, and it subsequently resamples the returned properties on
        the requested wavelength grid. */
        double getOpticalProperties(const Array& lambdav, const Array& thetav,
                                    Array& sigmaabsv, Array& sigmascav, Array& asymmparv,
                                    Table<2>& S11vv, Table<2>& S12vv, Table<2>& S33vv, Table<2>& S34vv) const override;

    /** This function must be implemented in each subclass to store the wavelengths and the
        corresponding tabulated properties in the array arguments, and to return the dust mass per
        hydrogen atom. The function must guarantee that all arrays have the same size. */
    virtual double getDustProperties(Array& lambdav, Array& kappaextv, Array& albedov, Array& asymmparv) const = 0;

    //======================== Other Functions =======================

public:
    /** This function returns the scattering mode supported by this material mix. For this class,
        it returns the HenyeyGreenstein mode. */
    ScatteringMode scatteringMode() const override;
};

////////////////////////////////////////////////////////////////////

#endif
