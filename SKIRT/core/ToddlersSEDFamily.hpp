/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TODDLERSSEDFAMILY_HPP
#define TODDLERSSEDFAMILY_HPP

#include "SEDFamily.hpp"
#include "StoredTable.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the ToddlersSEDFamily class represents the TODDLERS (Time evolution of Dust
    Diagnostics and Line Emission from Regions containing young Stars) family of star-forming
    region SEDs parameterized on age, metallicity, star formation efficiency, and cloud number
    density. This %SED template library is described in detail by Kapoor et al. in the manuscript
    submitted to MNRAS on May 26, 2023. The spectra were generated using the Cloudy master branch
    commit with SHA 69c3fa5871da3262341910e37c6ed2e5fb76dd3c.

    The SEDs in the TODDLERS library are tabulated over a wavelength range from \f$0.01\f$ to
    \f$3000~\mu\mathrm{m}\f$ and include both continuum and line emission. The low-resolution
    version of the library uses the default Cloudy output resolution of \f$R= \lambda/\Delta\lambda
    = 300\f$ to represent both the line and continuum components as calculated by Cloudy, resulting
    in 3787 wavelength points. The high-resolution version includes a large subset of the emission
    lines (as detailed in the paper), represented using high-resolution Gaussian line profiles with
    \f$R= \lambda/\Delta\lambda_\mathrm{G} = 5\times 10^4\f$, where \f$\Delta\lambda_\mathrm{G}\f$
    is the width of the \f$4 \sigma\f$ truncated Gaussian. These line profiles are added to the
    default \f$R=300\f$ resolution Cloudy continuum, resulting in a total of 9756 wavelength points
    with the overall spectral resolution structure shown below.

    \image html ToddlersSEDFamily.png

    The TODDLERS parameter space spans 90 ages from 0.1-30 Myr, 5 metallicity values from
    0.001-0.04, 7 star formation efficiencies from 0.01-0.15, and 9 cloud number densities from 
    10-2560 \f$\mathrm{cm}^{-3}\f$. The PAH-to-dust fraction is a maximum value scaled by the 
    neutral hydrogen abundance in component shells/clouds. Two maximum PAH fractions are available:
    High = 4.6%, Low = 1%.

    Finally, the SEDs are scaled by the mass of the imported star-forming region particle. The
    scaling assumes a cloud mass function obeying a power law with slope -1.8, i.e., \f$dN/dM
    \propto M^{-1.8}\f$ running from a cloud mass of \f$10^5\f$ to
    \f$10^{6.75}~\mathrm{M}_\odot\f$, sampled using 8 equi-logspaced cloud masses. Given the age
    range of the library, stellar particles younger than 30 Myr should be considered as
    star-forming region candidates, with their \em initial mass as the mass used for scaling.

    When imported from a text column file, the parameters for this %SED family must appear in the
    following order in the specified default units (unless these units are overridden by column
    header info): \f[ t\,(\mathrm{yr}) \quad Z\,(\mathrm{dimensionless}) \quad
    \mathrm{SFE}\,(\mathrm{dimensionless}) \quad n_\mathrm{cloud} (\mathrm{cm}^{-3}) \quad
    M\,(\mathrm{M}_\odot) \f] */
class ToddlersSEDFamily : public SEDFamily
{
    /** The enumeration type indicating the maximum PAH-to-dust fraction. */
    ENUM_DEF(PAHFraction, High, Low)
        ENUM_VAL(PAHFraction, High, "High PAH-to-dust fraction (4.6%)")
        ENUM_VAL(PAHFraction, Low, "Low PAH-to-dust fraction (1%)")
    ENUM_END()

    /** The enumeration type indicating the wavelength resolution. */
    ENUM_DEF(Resolution, Low, High)
        ENUM_VAL(Resolution, Low, "Low wavelength resolution (continuum and lines at R=300)")
        ENUM_VAL(Resolution, High, "High wavelength resolution (continuum at R=300 and lines at R=5e4)")
    ENUM_END()

    ITEM_CONCRETE(ToddlersSEDFamily, SEDFamily, "a Toddlers SED family for emission from star-forming regions")

        PROPERTY_ENUM(pahfraction, PAHFraction, "the maximum PAH-to-dust fraction")
        ATTRIBUTE_DEFAULT_VALUE(pahfraction, "High")

        PROPERTY_ENUM(resolution, Resolution, "the wavelength resolution")
        ATTRIBUTE_DEFAULT_VALUE(resolution, "Low")

    ITEM_END()

public:
    /** This constructor can be invoked programmatically by classes that use a hard-coded SED
        family (as opposed to selected through the ski file). Before the constructor returns, the
        newly created object is hooked up as a child to the specified parent in the simulation
        hierarchy (so it will automatically be deleted), and its setup() function has been called.
        */
    explicit ToddlersSEDFamily(SimulationItem* parent, PAHFraction pahfraction, Resolution resolution);

protected:
    /** This function opens the appropriate resource file (in SKIRT stored table format). */
    void setupSelfBefore() override;

    //====================== Other functions =====================

public:
    /** This function returns the number and type of parameters used by this particular %SED family
        as a list of SnapshotParameter objects. Each of these objects specifies unit information
        and a human-readable description for the parameter. */
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
