/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MAPPINGSSEDFAMILY_HPP
#define MAPPINGSSEDFAMILY_HPP

#include "SEDFamily.hpp"
#include "StoredTable.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the MappingsSEDFamily class represents the family of MAPPINGS III star-forming
    region template SEDs, parameterized on metallicity, compactness, ISM pressure and PDR covering
    factor, and scaled by star formation rate, as described by Groves et al. 2008 (ApJS,176,438).
    The model spectra have 1800 wavelength bins and cover 5 metallicities, 6 compactness
    parameters, 5 pressure values, and 2 covering factors (0 for pure HII and 1 for 100% PDR).

    The SEDs are tabulated over a wavelength range from 1 Angstrom to almost 1 meter
    with the spectral resolution shown in the figure below.

    \image html MappingsSEDFamily.png

    The data were downloaded from http://www.mpia-hd.mpg.de/~brent/starburst.html ->
    Cparam_models.save, however it seems that this page/download is no longer available. The
    contents of the IDL save file were converted to SKIRT stored table format for inclusion as a
    resource file. The stored table is opened during setup, and it is subsequently interpolated to
    the desired parameters and wavelength grid when needed.

    When imported from a text column file, the properties of the star-forming region represented by
    this %SED family must appear in the following order, and have the specified default units
    unless these units are overridden by column header info:

    \f[ \dot{M}\,(\mathrm{M}_\odot\,{\text{yr}}^{-1}) \quad Z\,(\mathrm{dimless}) \quad \log
    C\,(\mathrm{dimless}) \quad p\,(\mathrm{Pa}) \quad f_{\text{PDR}}\,(\mathrm{dimless}) \f]

    where \f$\dot{M}\f$ is the star formation rate, assumed to be constant over the past 10 Myr,
    \f$Z\f$ is the metallicity, \f$\log C\f$ is the logarithm of the compactness, \f$p\f$ is the
    ISM pressure, and \f$f_{\text{PDR}}\f$ is the PDR covering factor. */
class MappingsSEDFamily : public SEDFamily
{
    ITEM_CONCRETE(MappingsSEDFamily, SEDFamily, "a MAPPINGS III SED family for star-forming regions")
        ATTRIBUTE_TYPE_DISPLAYED_IF(MappingsSEDFamily, "Level2")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked programmatically by classes that use a hard-coded SED
        family (as opposed to selected through the ski file). Before the constructor returns, the
        newly created object is hooked up as a child to the specified parent in the simulation
        hierarchy (so it will automatically be deleted), and its setup() function has been called.
        */
    explicit MappingsSEDFamily(SimulationItem* parent);

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
    StoredTable<5> _table;
};

////////////////////////////////////////////////////////////////////

#endif
