/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RESOURCESED_HPP
#define RESOURCESED_HPP

#include "Array.hpp"
#include "ContSED.hpp"
#include "StoredTable.hpp"

////////////////////////////////////////////////////////////////////

/** ResourceSED is an abstract class for representing a spectral energy distribution loaded from a
    resource in SKIRT stored table format. The subclass must provide the name of the resource. */
class ResourceSED : public ContSED
{
    ITEM_ABSTRACT(ResourceSED, ContSED, "a spectral energy distribution loaded from a resource")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function opens the stored table resource tabulating the %SED and sets up the
        cumulative distribution that will be used to sample random wavelengths. */
    void setupSelfBefore() override;

    /** This function must be implemented in each subclass to return the name of the stored table
        resource tabulating the %SED. */
    virtual string resourceName() const = 0;

    //======================== Other Functions =======================

public:
    /** This function returns the intrinsic wavelength range of the %SED. For the current class,
        the range is retrieved from the underlying stored table. */
    Range intrinsicWavelengthRange() const override;

    /** This function returns the normalized specific luminosity \f$L_\lambda\f$ (i.e. radiative
        power per unit of wavelength) at the specified wavelength. */
    double specificLuminosity(double wavelength) const override;

    /** This function returns the normalized specific luminosity \f$L_\lambda\f$ (i.e. radiative
        power per unit of wavelength) at a number of wavelength points within the specified
        wavelength range. */
    void specificLuminosityArray(Array& lambdav, Array& pv, const Range& wavelengthRange) const override;

    /** This function returns the normalized integrated luminosity \f$L\f$ (i.e. radiative power)
        over the specified wavelength range. */
    double integratedLuminosity(const Range& wavelengthRange) const override;

    /** This function draws a random wavelength from the normalized spectral energy distribution.
        */
    double generateWavelength() const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    StoredTable<1> _table;
    Array _lambdav;
    Array _pv;
    Array _Pv;
    double _Ltot{0};
};

////////////////////////////////////////////////////////////////////

#endif
