/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MARASTONSEDFAMILY_HPP
#define MARASTONSEDFAMILY_HPP

#include "SEDFamily.hpp"
#include "StoredTable.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the MarastonSEDFamily class represents the family of Maraston SEDs for single
    stellar populations (SSPs), parameterized on metallicity and age and scaled by the initial mass
    of the SSP (Maraston 1998, MNRAS, 300, 872–892). We use both the Kroupa and Salpeter models,
    in each case for the red horizontal branch morphology (the models for the blue horizontal
    branch morphology cover only a fraction of the parameter space). In other words, the %SED
    family is available for two assumed initial mass functions (Kroupa and Salpeter).

    For metallicities [Z/H] = -1.35, -0.33, 0, 0.35, models are given for ages from 10^3 yr to 15
    Gyr. For metallicities [Z/H] = -2.25 and +0.67, only models from 1 to 15 Gyr are provided.
    These six metallicity values correspond to \f$Z \approx 0.00011247, 0.00089337, 0.0093547,
    0.02, 0.04477442, 0.09354703\f$.

    The SEDs are tabulated over a wavelength range from 0.009 \f$\mu\mathrm{m}\f$ to 160
    \f$\mu\mathrm{m}\f$ with the spectral resolution shown in the figure below.

    \image html MarastonSEDFamily.png

    The data were downloaded from
    http://www.icg.port.ac.uk/~maraston/Claudia's_Stellar_Population_Model.html and converted to
    SKIRT stored table format for inclusion as a resource file. The stored table is opened during
    setup, and it is subsequently interpolated to the desired parameters and wavelength grid when
    needed.

    When imported from a text column file, the parameters for this %SED family must appear in the
    following order in the specified default units (unless these units are overridden by column
    header info): \f[ M_\mathrm{init}\,(\mathrm{M}_\odot) \quad Z\,(\mathrm{dimensionless}) \quad
    t\,(\mathrm{yr}) \f] */
class MarastonSEDFamily : public SEDFamily
{
    /** The enumeration type indicating the assumed initial mass function (IMF). */
    ENUM_DEF(IMF, Kroupa, Salpeter)
        ENUM_VAL(IMF, Kroupa, "Kroupa IMF")
        ENUM_VAL(IMF, Salpeter, "Salpeter IMF")
    ENUM_END()

    ITEM_CONCRETE(MarastonSEDFamily, SEDFamily, "a Maraston SED family for single stellar populations")

        PROPERTY_ENUM(imf, IMF, "the assumed initial mass function")
        ATTRIBUTE_DEFAULT_VALUE(imf, "Kroupa")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked programmatically by classes that use a hard-coded SED
        family (as opposed to selected through the ski file). Before the constructor returns, the
        newly created object is hooked up as a child to the specified parent in the simulation
        hierarchy (so it will automatically be deleted), and its setup() function has been called.
        */
    explicit MarastonSEDFamily(SimulationItem* parent, IMF imf);

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
