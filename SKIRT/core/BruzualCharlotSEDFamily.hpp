/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BRUZUALCHARLOTSEDFAMILY_HPP
#define BRUZUALCHARLOTSEDFAMILY_HPP

#include "SEDFamily.hpp"
#include "StoredTable.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the BruzualCharlotSEDFamily class represents the family of Bruzual & Charlot
    2003 SEDs for single stellar populations (SSPs), parameterized on metallicity and age and
    scaled by the initial mass of the SSP (Bruzual & Charlot 2003, RAS 344, 1000-1026). We use the
    original 2003 Padova1994/Chabrier and Padova1994/Salpeter models, the two recommended models.
    In other words, the %SED family is available for two assumed initial mass functions (Chabrier
    and Salpeter). Each of those families can be loaded in two wavelength resolution versions, with
    respectively 1221 and 6900 wavelength points over a wavelength range from 0.009
    \f$\mu\mathrm{m}\f$ to 160 \f$\mu\mathrm{m}\f$. The spectral resolution for each version is
    shown in the figure below. The low-resolution version uses less resources and is sufficient for
    most purposes. However, if the model under study zooms in on narrow wavelength ranges in the
    optical regime, the high-resolution version may be preferable.

    \image html BruzualCharlotSEDFamily.png

    The data were downloaded from http://www.bruzual.org/~gbruzual/bc03/ and converted to SKIRT
    stored table format for inclusion as a resource file. The stored table is opened during setup,
    and it is subsequently interpolated to the desired parameters and wavelength grid when needed.

    When imported from a text column file, the parameters for this %SED family must appear in the
    following order in the specified default units (unless these units are overridden by column
    header info): \f[ M_\mathrm{init}\,(\mathrm{M}_\odot) \quad Z\,(\mathrm{dimensionless}) \quad
    t\,(\mathrm{yr}) \f] */
class BruzualCharlotSEDFamily : public SEDFamily
{
    /** The enumeration type indicating the assumed initial mass function (IMF). */
    ENUM_DEF(IMF, Chabrier, Salpeter)
        ENUM_VAL(IMF, Chabrier, "Chabrier IMF")
        ENUM_VAL(IMF, Salpeter, "Salpeter IMF")
    ENUM_END()

    /** The enumeration type indicating the wavelength resolution. */
    ENUM_DEF(Resolution, Low, High)
        ENUM_VAL(Resolution, Low, "Low wavelength resolution (1221 points)")
        ENUM_VAL(Resolution, High, "High wavelength resolution (6900 points)")
    ENUM_END()

    ITEM_CONCRETE(BruzualCharlotSEDFamily, SEDFamily, "a Bruzual-Charlot SED family for single stellar populations")

        PROPERTY_ENUM(imf, IMF, "the assumed initial mass function")
        ATTRIBUTE_DEFAULT_VALUE(imf, "Chabrier")

        PROPERTY_ENUM(resolution, Resolution, "the wavelength resolution")
        ATTRIBUTE_DEFAULT_VALUE(resolution, "Low")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked programmatically by classes that use a hard-coded SED
        family (as opposed to selected through the ski file). Before the constructor returns, the
        newly created object is hooked up as a child to the specified parent in the simulation
        hierarchy (so it will automatically be deleted), and its setup() function has been called.
        */
    explicit BruzualCharlotSEDFamily(SimulationItem* parent, IMF imf, Resolution resolution);

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
