/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FAMILYSED_HPP
#define FAMILYSED_HPP

#include "Array.hpp"
#include "ContSED.hpp"
class SEDFamily;

////////////////////////////////////////////////////////////////////

/** FamilySED is an abstract class for representing a spectral energy distribution obtained from an
    SED family for a set of parameters configured by the user. The subclass must provide a function
    to create the appropriate SEDFamily object and to return the configuration parameters. */
class FamilySED : public ContSED
{
    ITEM_ABSTRACT(FamilySED, ContSED, "a spectral energy distribution obtained from an SED family")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function asks the subclass to create the appropriate SEDFamily object and to return
        the configuration parameters. It then sets up the cumulative distribution that will be used
        to sample random wavelengths. */
    void setupSelfBefore() override;

    /** This function must be implemented in each subclass to return a newly created SEDFamily
        object (which is already hooked into the simulation item hierachy so it will be
        automatcially deleted) and to store the parameters for the specific %SED configured by the
        user in the specified array. The %SED will be normalized by this abstract base class, so the
        parameter values can be given for arbitary scaling. */
    virtual const SEDFamily* getFamilyAndParameters(Array& parameters) = 0;

    //======================== Other Functions =======================

public:
    /** This function returns the intrinsic wavelength range of the %SED. For the current class,
        the range is retrieved from the associated %SED family. */
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
    const SEDFamily* _family;
    Array _parameters;
    Array _lambdav;
    Array _pv;
    Array _Pv;
    double _Ltot{0};
};

////////////////////////////////////////////////////////////////////

#endif
