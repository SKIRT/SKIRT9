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
    includes binary stellar systems and that assumes a Chabrier IMF with an upper mass limit of 300
    solar masses (the "bin-imf_chab300" model).

    The SED templates are parametrized on metallicity (1e-5 - 0.04) and age (1 Myr - 100 Gyr), and
    scale with the initial mass of the SSP. The wavelength grid has a resolution of 1 Angstrom over
    the complete range (1 Angstrom - 10 micron). This corresponds to a spectral resolution as shown
    in the figure below.

    \image html BpassSEDFamily.png

    The data were downloaded from the BPASS web site at https://bpass.auckland.ac.nz/9.html and
    converted to SKIRT stored table format for inclusion in the <b>optional</b> BPASS resource
    pack. The stored table is opened during setup, and it is subsequently interpolated to the
    desired parameters and wavelength grid when needed.

    When imported from a text column file, the parameters for this %SED family must appear in the
    following order in the specified default units (unless these units are overridden by column
    header info): \f[ M_\mathrm{init}\,(\mathrm{M}_\odot) \quad Z\,(\mathrm{dimensionless}) \quad
    t\,(\mathrm{yr}) \f]

    <b>IMPORTANT NOTE TO THE USER</b>

    Because of its substantial size (450 MB compressed) the resource file required by this class is
    contained in a separate, <b>optional</b> resource pack, which can be installed by running the
    \c downloadResources.sh shell script in the SKIRT \c git directory or can be obtained from the
    download page on the SKIRT web site. */
class BpassSEDFamily : public SEDFamily
{
    ITEM_CONCRETE(BpassSEDFamily, SEDFamily, "a BPASS SED family for single stellar populations")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked programmatically by classes that use a hard-coded SED
        family (as opposed to selected through the ski file). Before the constructor returns, the
        newly created object is hooked up as a child to the specified parent in the simulation
        hierarchy (so it will automatically be deleted), and its setup() function has been called.
        */
    explicit BpassSEDFamily(SimulationItem* parent);

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
