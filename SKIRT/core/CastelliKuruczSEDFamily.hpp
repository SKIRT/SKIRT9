/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CASTELLIKURUCZSEDFAMILY_HPP
#define CASTELLIKURUCZSEDFAMILY_HPP

#include "SEDFamily.hpp"
#include "StoredTable.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the CastelliKuruczSEDFamily class represents the Castelli-Kurucz atlas of
    stellar atmosphere models for a wide range of metallicities, effective temperatures, and
    gravities. These models were computed by Fiorella Castelli (Castelli & Kurucz 2003, IAU
    Symposium 210, Modelling of Stellar Atmospheres, Uppsala, Sweden, eds. N.E. Piskunov, W.W.
    Weiss. and D.F. Gray, 2003, ASP-S210) and include improvements upon those previously provided
    by by Kurucz (1990).

    The luminosity values given by the models are normalized for a stellar surface equal to the
    unit sphere, so these values must be multipled by the stellar surface, i.e. \f$4\pi\,R^2\f$.

    The parameter space is discretized on 8 metallicities, 76 temperatures, and 11 gravity values
    for a total of 6688 possible models. From these possibilities, 3808 models are actually
    available and 2880 have not been calculated. The parameter ranges for available models are
    summarized in the table below, with in all cases \f$ 6.325 \times 10^{-5} \le Z \le 6.325
    \times 10^{-2} \f$.

    | Effective temperature (K) | | Gravity (\f$\log_{10} (g/\mathrm{cm}\,\mathrm{s}^{-2})\f$) |
    |----|----|----|
    |  3500 <= Teff <=  6000    | --> |    0.0 <= logg <= 5.0   |
    |  6000 <  Teff <=  7500    | --> |    0.5 <= logg <= 5.0   |
    |  7500 <  Teff <=  8250    | --> |    1.0 <= logg <= 5.0   |
    |  8250 <  Teff <=  9000    | --> |    1.5 <= logg <= 5.0   |
    |  9000 <  Teff <= 11750    | --> |    2.0 <= logg <= 5.0   |
    | 11750 <  Teff <= 19000    | --> |    2.5 <= logg <= 5.0   |
    | 19000 <  Teff <= 26000    | --> |    3.0 <= logg <= 5.0   |
    | 26000 <  Teff <= 31000    | --> |    3.5 <= logg <= 5.0   |
    | 31000 <  Teff <= 39000    | --> |    4.0 <= logg <= 5.0   |
    | 39000 <  Teff <= 49000    | --> |    4.5 <= logg <= 5.0   |
    | 49000 <  Teff <= 50000    | --> |           logg  = 5.0   |

    The following table lists the properties of solar metallicity stars (\f$Z=0.02\f$) of different
    spectral types and luminosity classes as suggested by Castelli and Kurucz based on other sources.

    | Stellar type | Temperature (K) | Gravity (\f$\log_{10} (g/\mathrm{cm}\,\mathrm{s}^{-2})\f$) |
    |----|----|----|
    |  O3V    |  44852 |    +3.92
    |  O5V    |  40862 |    +3.92
    |  O5.5V  |  39865 |    +3.92
    |  O6V    |  38867 |    +3.92
    |  O6.5V  |  37870 |    +3.92
    |  O7V    |  36872 |    +3.92
    |  O7.5V  |  35874 |    +3.92
    |  O8V    |  34877 |    +3.92
    |  O8.5   |  33879 |    +3.92
    |  O9V    |  32882 |    +3.92
    |  O9.5   |  31884 |    +3.92
    |  B0V    |  30000 |    +3.90
    |  B1V    |  25400 |    +3.90
    |  B3V    |  18700 |    +3.94
    |  B5V    |  15400 |    +4.04
    |  B8V    |  11900 |    +4.04
    |  A0V    |   9520 |    +4.14
    |  A1V    |   9230 |    +4.10
    |  A3V    |   8270 |    +4.20
    |  A5V    |   8200 |    +4.29
    |  F0V    |   7200 |    +4.34
    |  F2V    |   6890 |    +4.34
    |  F5V    |   6440 |    +4.34
    |  F8V    |   6200 |    +4.40
    |  G0V    |   6030 |    +4.39
    |  G2V    |   5860 |    +4.40
    |  G8V    |   5570 |    +4.50
    |  K0V    |   5250 |    +4.49
    |  K2V    |   4780 |     +4.5
    |  K4V    |   4560 |     +4.5
    |  K5V    |   4350 |    +4.54
    |  K7V    |   4060 |     +4.5
    |  M0V    |   3850 |    +4.59
    |  M2V    |   3580 |    +4.64
    |  M6V    |   3050 |    +5.00
    |  B0III  |  29000 |    +3.34
    |  B5III  |  15000 |    +3.49
    |  G0III  |   5850 |    +2.94
    |  G5III  |   5150 |    +2.54
    |  K0III  |   4750 |    +2.14
    |  K5III  |   3950 |    +1.74
    |  M0III  |   3800 |    +1.34
    |  BOI    |  26000 |    +2.84
    |  B5I    |  13600 |    +2.44
    |  AOI    |   9730 |    +2.14
    |  A5I    |   8510 |    +2.04
    |  F0I    |   7700 |    +1.74
    |  F5I    |   6900 |    +1.44
    |  G0I    |   5550 |    +1.34
    |  G5I    |   4850 |    +1.14
    |  K0I    |   4420 |    +0.94
    |  K5I    |   3850 |    +0.00
    |  M0I    |   3650 |    -0.10
    |  M2I    |   3600 |    -0.10

    The SEDs are tabulated over a wavelength range from 0.009 \f$\mu\mathrm{m}\f$ to 160
    \f$\mu\mathrm{m}\f$ with the spectral resolution shown in the figure below.

    \image html CastelliKuruczSEDFamily.png

    The data for the models were downloaded from
    http://www.stsci.edu/hst/observatory/crds/castelli_kurucz_atlas.html and converted to SKIRT
    stored table format for inclusion as a resource file. The stored table is opened during setup,
    and it is subsequently interpolated to the desired parameters and wavelength grid when needed.

    When imported from a text column file, the properties of the star-forming region represented by
    this %SED family must appear in the following order, and have the specified default units
    unless these units are overridden by column header info:

    \f[ R\,(\mathrm{km}) \quad Z\,(\mathrm{dimless}) \quad T_\mathrm{eff}\,(\mathrm{K})
    \quad g\,(\mathrm{m}\,\mathrm{s}^{-2}) \f]

    where \f$R\f$ is the stellar radius, \f$T_\mathrm{eff}\f$ is the effective stellar surface
    temperature, \f$Z\f$ is the stellar metallicity, and \f$g\f$ is the gravity at the stellar
    surface. */
class CastelliKuruczSEDFamily : public SEDFamily
{
    ITEM_CONCRETE(CastelliKuruczSEDFamily, SEDFamily, "a Castelli-Kurucz SED family for stellar atmospheres")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked programmatically by classes that use a hard-coded SED
        family (as opposed to selected through the ski file). Before the constructor returns, the
        newly created object is hooked up as a child to the specified parent in the simulation
        hierarchy (so it will automatically be deleted), and its setup() function has been called.
        */
    explicit CastelliKuruczSEDFamily(SimulationItem* parent);

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
    StoredTable<4> _table;
};

////////////////////////////////////////////////////////////////////

#endif
