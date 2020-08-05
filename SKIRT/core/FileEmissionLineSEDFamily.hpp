/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILEEMISSIONLINESEDFAMILY_HPP
#define FILEEMISSIONLINESEDFAMILY_HPP

#include "SEDFamily.hpp"
#include "StoredTable.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the FileEmissionLineSEDFamily can be used to import user-provided spectra into a SKIRT simulation.
    This is for instance useful when including line emission from the ISM itself in an AdaptiveMeshSource, CellSource,
    ParticleSource...

    The spectra should be provided in SKIRT stored table format as a 2-dimensional table with the wavelength as the
    first axis.

    The implementation has on purpose been kept very basic to allow for maximal flexibility: the spectra are
    parametrised using a single (positive) dimensionless index parameter that identifies each spectrum. Apart from
    the positivity requirement, there are no restrictions on the values of the index. Although originally conceived
    as an integer counter number, the use of a stored table by the implementation also allows for floating point
    indices. Care has to be taken that the indices referenced by the Source do in fact exist; the stored table will
    (wrongly) interpolate the spectra for non-existing indices.

    To allow for cells without emission, negative indices are ignored. Both specificLuminosity() and cdf() return 0
    whenever a negative input index is specified. */
class FileEmissionLineSEDFamily : public SEDFamily
{
    ITEM_CONCRETE(FileEmissionLineSEDFamily, SEDFamily, "a user-provided list of emission lines that make up an SED")
        PROPERTY_STRING(filename, "the name of the stored table file defining the emission lines")
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
