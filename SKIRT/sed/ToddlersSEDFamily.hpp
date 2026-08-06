/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TODDLERSSEDFAMILY_HPP
#define TODDLERSSEDFAMILY_HPP

#include "SEDFamily.hpp"
#include "StoredTable.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the ToddlersSEDFamily class represents the family of star-forming region
    templates from the TODDLERS model suite (Time evolution of Observables including Dust
    Diagnostics and Line Emission from Regions containing young Stars).

    The TODDLERS model calculates the spherical evolution of a gas cloud around a young stellar
    cluster, accounting for stellar feedback processes such as stellar winds, supernovae, radiation
    pressure, and gravitational forces; see Kapoor et al. 2023 (MNRAS, 526, 3871) and Kapoor et al.
    2024 (A&A, 692A, 79). The spectra were generated using the third-party Cloudy code based on the
    results of the evolution model.

    <b>%Configuration Options</b>

    This class offers several user-configured options as described below.

    \em sedMode determines how the %SED is calculated and scaled:
      - SFRNormalized: SEDs pre-integrated over time and cloud mass spectrum, directly scaled by SFR
      - Cloud: SEDs for individual star-forming clouds with explicit time evolution

    \em stellarTemplate determines the stellar population model, IMF, and stellar evolution:
      - SB99Kroupa100Sin: Starburst99 models with Kroupa IMF (0.1-100 \f$\mathrm{M}_\odot\f$) and
        single star evolution
      - BPASSChab100Bin: BPASS models with Chabrier IMF (0.1-100 \f$\mathrm{M}_\odot\f$) and binary
        star evolution
      - BPASSChab300Bin: BPASS models with Chabrier IMF (0.1-300 \f$\mathrm{M}_\odot\f$) and binary
        star evolution

    \em includeDust determines dust processing in the SEDs:

    - true:
        - For low spectral resolution, uses the total emission reported by Cloudy, which includes
          attenuated stellar, nebular, and dust continuum components.

        - For high spectral resolution, uses dust-attenuated emission lines (emergent luminosities)
          added to the total emission after removing the low-resolution lines. The high-resolution
          mode includes a limited set of approximately 140 emission lines tracked by TODDLERS, so the
          replacement is not one-to-one. Emergent luminosity values are calculated using escape
          probabilities for diffuse radiation in Cloudy. Roughly, these drop off as \f$E_2(\tau)\f$,
          where \f$E_2\f$ is the exponential integral and \f$\tau\f$ is the optical depth.
          \f$E_2(\tau)\f$ drops off faster than \f$\exp(-\tau)\f$. This is discussed in section 3.1.1
          of Rollig et al. 2007 (A&A, 467, 187).

    - false:
        - For low spectral resolution, uses only the incident stellar radiation field (stellar
          continuum) without any gas or dust processing.
        - For high spectral resolution, uses intrinsic emission line luminosities (without foreground
          attenuation) added to the incident stellar continuum.

    \em resolution determines the spectral resolution of the SEDs:

    - Low: the entire spectrum (continuum and lines) is at \f$R=300\f$ resolution

    - High: low-resolution continuum with selected emission lines represented as high-resolution
    Gaussian profiles (\f$R = \lambda/\Delta\lambda_G = 5\times 10^4\f$) sampled using 37 points
    per line. This approach offers better line-to-continuum contrast while maintaining
    computational efficiency.

    \em sfrPeriod sets the time period over which SFR is averaged and integrated (only used in
    SFRNormalized mode):

    - Period10Myr: SFR integrated over 10 Myr
    - Period30Myr: SFR integrated over 30 Myr

    <b>Wavelength Coverage</b>

    The SEDs in the TODDLERS library are tabulated over a wavelength range from \f$0.01\f$ to
    \f$3000~\mu\mathrm{m}\f$ (UV through millimeter) and include stellar, nebular, and dust
    continuum emission along with numerous emission lines from H II, PDR, and molecular gas phases
    from Cloudy spectral synthesis calculations. For more details, see the discussion on dust
    processing and spectral resolution in the previous section. The figure below shows the spectral
    resolution for a high-resolution model.

    \image html ToddlersSEDFamily.png

    <em>Emission Line Artifacts</em>

    For the high-resolution Cloud-mode models with includeDust=true, the process of subtracting
    low-resolution line emission from the total low-resolution spectrum before adding
    high-resolution line profiles can occasionally produce numerical artifacts. These artifacts
    appear as sharp drops in flux sometimes up to 5 orders of magnitude below the continuum level.
    This has been observed particularly around mid-IR lines in the 10-30 \f$\mu\mathrm{m}\f$
    region. The artifacts primarily affect lower mass clouds (\f$\le 10^{5.5}~\mathrm{M}_\odot\f$)
    and occur at the transition points between the continuum and line profiles. The artifacts are
    expected to have a negligible impact on the integrated line luminosities in actual galaxy
    simulations using an ensemble of clouds.

    <b>Parameter space</b>

   The SEDs in the TODDLERS library are tabulated for combinations of the following parameters.

    1. Evolution time (Cloud mode only): the time since the start of evolution
       - Time-dependent evolution from 0.1 to 30 Myr (90 values)
    2. Metallicity:
       - Starburst99: range from Z=0.001 to Z=0.04 (5 values)
       - BPASS: range from Z=0.001 to Z=0.04 (11 values)
    3. Star formation efficiency (SFE): the fraction of cloud mass converted to stars
       - Starburst99: range from 0.01 to 0.15 (7 values)
       - BPASS: range from 0.01 to 0.1 (5 values)
    4. Cloud number density: the initial density of the star-forming cloud
       - Starburst99: range from 10 to 2560 \f$\mathrm{cm}^{-3}\f$ (9 values)
       - BPASS: range from 40 to 640 \f$\mathrm{cm}^{-3}\f$ (5 values)
    5. Cloud mass (Cloud mode only): the mass of the star-forming cloud
       - range from \f$10^5\f$ to \f$10^{6.75}~\mathrm{M}_\odot\f$ (8 values)

    When using SFRNormalized mode, the parameters must appear in the following order, with the
    specified default units unless overridden by column header info:

    \f[ Z\,(\mathrm{dimensionless}) \quad \mathrm{SFE}\,(\mathrm{dimensionless}) \quad
    n_{\text{cl}}\,(\mathrm{cm}^{-3}) \quad \dot{M}_*\,(\mathrm{M}_\odot\,\mathrm{yr}^{-1}) \f]

    where \f$Z\f$ is the metallicity, \f$\mathrm{SFE}\f$ is the star formation efficiency,
    \f$n_{\text{cl}}\f$ is the cloud number density, and \f$\dot{M}_*\f$ is the star formation
    rate.

    When using Cloud mode, the parameters must appear in the following order, with the specified
    default units unless overridden by column header info:

    \f[ t\,(\mathrm{Myr}) \quad Z\,(\mathrm{dimensionless}) \quad
    \mathrm{SFE}\,(\mathrm{dimensionless}) \quad n_{\text{cl}}\,(\mathrm{cm}^{-3}) \quad
    M_{\text{cl}}\,(\mathrm{M}_\odot) \quad \mathrm{scaling}\,(\mathrm{dimensionless}) \f]

    where \f$t\f$ is the evolution time since the start of star formation, \f$Z\f$ is the
    metallicity, \f$\mathrm{SFE}\f$ is the star formation efficiency, \f$n_{\text{cl}}\f$ is the
    cloud number density, \f$M_{\text{cl}}\f$ is the cloud mass, and \f$\mathrm{scaling}\f$ is an
    arbitrary scaling factor (typically 1).

    <b>Recollapse</b>

    If stellar feedback is insufficient to overcome gravity, the evolution of a star-forming cloud
    includes a recollapse phase, triggering a subsequent generation of star formation.
    Consequently, at a given point in time, multiple stellar populations of different ages may be
    present.

    In SFRNormalized mode, a constant star formation history is assumed over the past 10 or 30 Myr
    period (as determined by the \em sfrPeriod parameter). The model properly accounts for
    recollapse: any recollapse contribution is pre-integrated over the time evolution (10 or 30
    Myr) and cloud mass spectrum (\f$10^5\f$ to \f$10^{6.75}~\mathrm{M}_\odot\f$) with a power-law
    distribution (\f$dN/dM \propto M^{-1.8}\f$).

    In Cloud mode, users must apply the appropriate scaling to ensure that the star formation rate
    follows the input model when recollapse occurs. Assuming SKIRT is used to post-process a
    galaxy snapshot from a hydrodynamical simulation, this can be achieved with the following steps:

    1. Divide the young stellar population particles into temporal bins based on their age.
    2. Calculate recollapse contributions for each temporal bin.
    3. Apply a uniform scaling factor to all particles within a time bin based on the total
       recollapse contribution in that bin, thereby maintaining the simulated galaxy's total
       star formation rate.

    A Python implementation of this algorithm can be found at https://github.com/anandutsavkapoor/reSample

    In other words, when a recollapse event occurs, stellar mass is effectively added to one part
    of the galaxy while being removed from other parts based solely on the time bin in which the
    recollapse occurs, without any consideration of spatial proximity or physical connection
    between these regions. Although this approach preserves the galaxy's total star formation rate,
    users should consider this limitation when interpreting results, particularly for analyses that
    are sensitive to local star formation or that examine spatially resolved properties of
    simulated galaxies. The SFRNormalized mode circumvents this issue at the expense of explicit
    time evolution.

    <b>Resource packs</b>

    The tabulated data used by the various TODDLERS models is contained in SKIRT resource packs
    that need to be downloaded separately. Because of the large size of the Cloud-mode tables,
    the data are distributed over multiple resource packs. It is recommended to download and
    install only the packs for the required models:

    - SKIRT9_Resources_TODDLERS (400 MB installed):
            all SFRNormalized-mode models
    - SKIRT9_Resources_TODDLERS_Cloud_SB99_kroupa100 (50 GB installed):
            Cloud-mode models with Starburst99 and Kroupa IMF (0.1-100)
    - SKIRT9_Resources_TODDLERS_Cloud_BPASS_chab100 (43 GB installed):
            Cloud-mode models with BPASS and Chabrier IMF (0.1-100)
    - SKIRT9_Resources_TODDLERS_Cloud_BPASS_chab300 (43 GB installed):
            Cloud-mode models with BPASS and Chabrier IMF (0.1-300)

    */
class ToddlersSEDFamily : public SEDFamily
{
    /** The enumeration type indicating the SED calculation mode */
    ENUM_DEF(SedMode, SFRNormalized, Cloud)
        ENUM_VAL(SedMode, SFRNormalized, "SEDs normalized by star formation rate")
        ENUM_VAL(SedMode, Cloud, "Individual cloud SEDs with time evolution")
    ENUM_END()

    /** The enumeration type indicating the stellar template to use */
    ENUM_DEF(StellarTemplate, SB99Kroupa100Sin, BPASSChab100Bin, BPASSChab300Bin)
        ENUM_VAL(StellarTemplate, SB99Kroupa100Sin,
                 "Starburst99 with Kroupa IMF (0.1-100 Msun) and single star evolution")
        ENUM_VAL(StellarTemplate, BPASSChab100Bin, "BPASS with Chabrier IMF (0.1-100 Msun) and binary star evolution")
        ENUM_VAL(StellarTemplate, BPASSChab300Bin, "BPASS with Chabrier IMF (0.1-300 Msun) and binary star evolution")
    ENUM_END()

    /** The enumeration type indicating the wavelength resolution */
    ENUM_DEF(Resolution, Low, High)
        ENUM_VAL(Resolution, Low, "Low wavelength resolution (continuum and lines at R=300)")
        ENUM_VAL(Resolution, High, "High wavelength resolution (continuum at R=300 and lines at R=5e4)")
    ENUM_END()

    /** The enumeration type indicating the SFR integration period */
    ENUM_DEF(SFRPeriod, Period10Myr, Period30Myr)
        ENUM_VAL(SFRPeriod, Period10Myr, "SFR integrated over 10 Myr")
        ENUM_VAL(SFRPeriod, Period30Myr, "SFR integrated over 30 Myr")
    ENUM_END()

    ITEM_CONCRETE(ToddlersSEDFamily, SEDFamily, "a TODDLERS SED family for emission from star-forming regions")

        PROPERTY_ENUM(sedMode, SedMode, "SED calculation mode")
        ATTRIBUTE_DEFAULT_VALUE(sedMode, "SFRNormalized")

        PROPERTY_ENUM(stellarTemplate, StellarTemplate, "the stellar template, IMF, and evolution model to use")
        ATTRIBUTE_DEFAULT_VALUE(stellarTemplate, "SB99Kroupa100Sin")

        PROPERTY_BOOL(includeDust, "include dust processing in the SED models")
        ATTRIBUTE_DEFAULT_VALUE(includeDust, "true")

        PROPERTY_ENUM(resolution, Resolution, "the wavelength resolution")
        ATTRIBUTE_DEFAULT_VALUE(resolution, "Low")

        PROPERTY_ENUM(sfrPeriod, SFRPeriod, "the SFR integration time period")
        ATTRIBUTE_DEFAULT_VALUE(sfrPeriod, "Period10Myr")
        ATTRIBUTE_RELEVANT_IF(sfrPeriod, "sedModeSFRNormalized")

    ITEM_END()

public:
    /** This constructor can be invoked programmatically by classes that use a hard-coded SED
        family. The newly created object is hooked up as a child to the specified parent in the
        simulation hierarchy, and its setup() function has been called. */
    explicit ToddlersSEDFamily(SimulationItem* parent, SedMode sedMode, StellarTemplate stellarTemplate,
                               bool includeDust, Resolution resolution, SFRPeriod sfrPeriod = SFRPeriod::Period10Myr);

protected:
    /** This function opens the appropriate resource file (in SKIRT stored table format). */
    void setupSelfBefore() override;

public:
    /** This function returns the number and type of parameters used by this particular %SED family
        as a list of SnapshotParameter objects. Each object specifies unit information and a
        human-readable description for the parameter. */
    vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns the intrinsic wavelength range of the %SED family from the stored
        table. */
    Range intrinsicWavelengthRange() const override;

    /** This function returns the specific luminosity \f$L_\lambda\f$ (radiative power per unit of
        wavelength) for the %SED with the specified parameters at the specified wavelength. */
    double specificLuminosity(double wavelength, const Array& parameters) const override;

    /** This function constructs the normalized probability density function (pdf) and cumulative
        distribution function (cdf) for the %SED with the specified parameters. */
    double cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
               const Array& parameters) const override;

    //======================== Data Members =======================

private:
    // Only one of these tables will be used, depending on the sedMode
    StoredTable<6> _cloudTable;          // 6D table: lambda, time, Z, SFE, n_cl, M_cl
    StoredTable<4> _sfrNormalizedTable;  // 4D table: lambda, Z, SFE, n_cl
};

//////////////////////////////////////////////////////////////////////

#endif
