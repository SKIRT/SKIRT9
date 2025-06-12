/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BPASSSEDFAMILY_HPP
#define BPASSSEDFAMILY_HPP

#include "SEDFamily.hpp"
#include "StoredTable.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the BpassSEDFamily class represents a BPASS family of single stellar population
    (SSP) %SED templates (Eldridge, Stanway et al, 2017, PASA 34, 58; Stanway and Eldridge, 2018,
    MNRAS, 479, 75). We use the BPASS data release version 2.2.1 (July 2018) of the model that
    includes binary stellar systems and that assumes a Chabrier IMF.

    The SED templates are parametrized on metallicity (1e-5 - 0.04) and age (1 Myr - 100 Gyr), and
    scale with the initial mass of the SSP.

    This class offers two user-configured options:

    \em imf determines the IMF, more specifically the upper limit of the stellar mass range:
      - Chabrier100: Chabrier IMF with stellar masses from 0.1 to 100 \f$\mathrm{M}_\odot\f$
      - Chabrier300: Chabrier IMF with stellar masses from 0.1 to 300 \f$\mathrm{M}_\odot\f$

    \em resolution determines the spectral resolution and the wavelength coverage:
      - Original: the wavelength grid has a resolution of 1 Angstrom over the wavelength range
        between 1 Angstrom and 10 micron. This is the original resolution and wavelength coverage
        of the BPASS library, downloaded from the BPASS web site at
        https://bpass.auckland.ac.nz/9.html.
      - Downsampled: the wavelength coverage is downsampled. Between 1 Angstrom and 0.6
        micron, the original resolution of 1 Angstrom is retained. Between 0.6 and 10 micron the
        resolution is lowered to \f$R=711\f$ (2000 grid points). Finally, a flux point is added at 160
        micron to ensure that the Rayleigh-Jeans tail of the SED is properly covered. This one
        additional point is using the slope between 8 and 10 micron.

    When imported from a text column file, the parameters for this %SED family must appear in the
    following order in the specified default units (unless these units are overridden by column
    header info): \f[ M_\mathrm{init}\,(\mathrm{M}_\odot) \quad Z\,(\mathrm{dimensionless}) \quad
    t\,(\mathrm{yr}) \f]

    <b>IMPORTANT NOTE TO THE USER</b>

    Because of their substantial size, the resource files required by this class are
    contained in a separate, <b>optional</b> resource pack, which can be installed by running the
    \c downloadResources.sh shell script in the SKIRT \c git directory or can be obtained from the
    download page on the SKIRT web site. */

class BpassSEDFamily : public SEDFamily
{
    /** The enumeration type indicating the IMF to use */
    ENUM_DEF(IMF, Chabrier100, Chabrier300)
        ENUM_VAL(IMF, Chabrier100, "Chabrier IMF (0.1-100 Msun)")
        ENUM_VAL(IMF, Chabrier300, "Chabrier IMF (0.1-300 Msun)")
    ENUM_END()

    /** The enumeration type indicating the wavelength resolution and spectral range */
    ENUM_DEF(Resolution, Original, Downsampled)
        ENUM_VAL(Resolution, Original, "Original wavelength resolution and range")
        ENUM_VAL(Resolution, Downsampled, "Downsampled wavelength resolution and extended wavelength range")
    ENUM_END()

    ITEM_CONCRETE(BpassSEDFamily, SEDFamily, "a BPASS SED family for single stellar populations")

        PROPERTY_ENUM(imf, IMF, "IMF")
        ATTRIBUTE_DEFAULT_VALUE(imf, "Chabrier300")

        PROPERTY_ENUM(resolution, Resolution, "the wavelength resolution")
        ATTRIBUTE_DEFAULT_VALUE(resolution, "Original")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked programmatically by classes that use a hard-coded SED
        family (as opposed to selected through the ski file). Before the constructor returns, the
        newly created object is hooked up as a child to the specified parent in the simulation
        hierarchy (so it will automatically be deleted), and its setup() function has been called.
        */
    explicit BpassSEDFamily(SimulationItem* parent, IMF imf, Resolution resolution);

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
