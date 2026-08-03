/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BLACKBODYSEDFAMILY_HPP
#define BLACKBODYSEDFAMILY_HPP

#include "SEDFamily.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the BlackBodySEDFamily class represents a family of black body %SEDs, i.e.
    emission spectra of perfect absorbers in thermal equilibrium. The %SED family is parameterized
    on the temperature \f$T\f$ of the black body (defining the form of the spectrum) and the radius
    \f$R\f$ of the black body (determining the absolute luminosity scale). As a function of these
    two parameters, the specific luminosity is given by \f[ L_\lambda(\lambda,R,T) =
    4\pi\,R^2\;\pi\;B_\lambda(\lambda,T) \f] where \f$B_\lambda(\lambda,T)\f$ is the Planck
    function.

    Whenever a tabular form of a black body %SED is requested, this class uses a spectral
    resolution of \f$R\triangleq\lambda/\Delta\lambda\geq 1000\f$ (see PlanckFunction::cdf()).

    When imported from a text column file, the parameters for this %SED family must appear in the
    following order in the specified default units (unless these units are overridden by column
    header info): \f[ R\,(\mathrm{km}) \quad T\,(\mathrm{K}) \f] */
class BlackBodySEDFamily : public SEDFamily
{
    ITEM_CONCRETE(BlackBodySEDFamily, SEDFamily, "a black body SED family")
    ITEM_END()

    //====================== Other functions =====================

public:
    /** This function returns the number and type of parameters used by this particular %SED family
        as a list of SnapshotParameter objects. Each of these objects specifies unit information
        and a human-readable descripton for the parameter. */
    vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns the intrinsic wavelength range of the %SED family. For the
        BlackBodySEDFamily, the intrinsic range is unlimited, so this function returns a range
        including all representable positive floating point values. */
    Range intrinsicWavelengthRange() const override;

    /** This function returns the specific luminosity \f$L_\lambda\f$ (i.e. radiative power per
        unit of wavelength) for the %SED with the specified parameters at the specified wavelength,
        or zero if the wavelength is outside of the %SED's intrinsic wavelength range. The number
        and type of parameters must match the information returned by the parameterInfo() function;
        if not the behavior is undefined. */
    double specificLuminosity(double wavelength, const Array& parameters) const override;

    /** This function constructs both the normalized probability density function (pdf) and the
        corresponding normalized cumulative distribution function (cdf) for the %SED with the
        specified parameters over the specified wavelength range. The function returns the
        normalization factor. The number and type of parameters must match the information returned
        by the parameterInfo() function; if not the behavior is undefined. */
    double cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
               const Array& parameters) const override;
};

////////////////////////////////////////////////////////////////////

#endif
