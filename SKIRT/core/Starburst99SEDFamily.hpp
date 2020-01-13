/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STARBURST99SEDFAMILY_HPP
#define STARBURST99SEDFAMILY_HPP

#include "SEDFamily.hpp"
#include "StoredTable.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the Starburst99SEDFamily class represents a family of Starburst99 SEDs for
    single stellar populations (Leitherer et al. 1999, ApJS, 123, 3), assuming the Kroupa initial
    mass function (Kroupa 2001, MNRAS, 322, 231), parameterized on metallicity and age, and scaled
    by the initial mass of the SSP. The library data was prepared and bundled into a FITS file by
    Patrik Jonsson for use by the \c Sunrise code (<tt>2013ascl.soft03030J</tt>).

    The SEDs are tabulated over a wavelength range from 0.009 \f$\mu\mathrm{m}\f$ to 160
    \f$\mu\mathrm{m}\f$ with the spectral resolution shown in the figure below.

    \image html Starburst99SEDFamily.png

    The data were downloaded from
    <tt>https://bitbucket.org/lutorm/sunrise/downloads/Patrik-imfKroupa-Zmulti-ml.fits.gz</tt> and
    converted to SKIRT stored table format for inclusion as a resource file. The stored table is
    opened during setup, and it is subsequently interpolated to the desired parameters and
    wavelength grid when needed.

    When imported from a text column file, the parameters for this %SED family must appear in the
    following order in the specified default units (unless these units are overridden by column
    header info): \f[ M_\mathrm{init}\,(\mathrm{M}_\odot) \quad Z\,(\mathrm{dimensionless}) \quad
    t\,(\mathrm{yr}) \f] */
class Starburst99SEDFamily : public SEDFamily
{
    ITEM_CONCRETE(Starburst99SEDFamily, SEDFamily, "a Starburst99 SED family for single stellar populations")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked programmatically by classes that use a hard-coded SED
        family (as opposed to selected through the ski file). Before the constructor returns, the
        newly created object is hooked up as a child to the specified parent in the simulation
        hierarchy (so it will automatically be deleted), and its setup() function has been called.
        */
    explicit Starburst99SEDFamily(SimulationItem* parent);

protected:
    /** This function opens the appropriate resource file (in SKIRT stored table format). */
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
    StoredTable<3> _table;
};

////////////////////////////////////////////////////////////////////

#endif
