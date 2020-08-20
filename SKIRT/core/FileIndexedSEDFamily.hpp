/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILEINDEXEDSEDFAMILY_HPP
#define FILEINDEXEDSEDFAMILY_HPP

#include "SEDFamily.hpp"
#include "StoredTable.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the FileIndexedSEDFamily can be used to specify user-provided spectra for an
    imported source such as a CellSource or a ParticleSource. This can be useful, for example, to
    include emission from the interstellar medium in a simulation. The spectra should be provided
    in SKIRT stored table format as a 2-dimensional table with the wavelength as the first axis.

    The implementation of this class has on purpose been kept very basic to allow for maximum
    flexibility: the spectra are parametrised using a single (positive) dimensionless index
    parameter that identifies each spectrum. Apart from the positivity requirement, there are no
    restrictions on the values of the index. Although originally conceived as an integer counter,
    the use of a stored table by the implementation also allows for floating point indices. Care
    has to be taken that the indices referenced by the imported source elements do in fact exist in
    the stored table. If not, the implementation will interpolate from existing index values and
    this interpolation will most likely be wrong because spectra with consecutive indices do not
    necessarily have a physical relationship.

    To allow for source elements without emission, a negative index is interpreted as "no
    emission". */
class FileIndexedSEDFamily : public SEDFamily
{
    ITEM_CONCRETE(FileIndexedSEDFamily, SEDFamily, "a user-provided, indexed SED family")
        PROPERTY_STRING(filename, "the name of the stored table file listing the SEDs")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function opens the user-specified resource file (in SKIRT stored table format). */
    void setupSelfBefore() override;

    //====================== Other functions =====================

public:
    /** This function returns the number and type of parameters used by this particular %SED family
        as a list of SnapshotParameter objects. Each of these objects specifies unit information
        and a human-readable descripton for the parameter. */
    vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns the intrinsic wavelength range of the %SED family. It retrieves this
        range from the underlying stored table. */
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

    //====================== Data members =====================

private:
    StoredTable<2> _table;
};

////////////////////////////////////////////////////////////////////

#endif
