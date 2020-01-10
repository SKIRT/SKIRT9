/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BROADBAND_HPP
#define BROADBAND_HPP

#include "Band.hpp"
#include "StoredTable.hpp"

////////////////////////////////////////////////////////////////////

/** A BroadBand object represents a standard wavelength band with a transmission curve that is
    loaded from a resource file provided with SKIRT. Examples include the standard Johnson filters
    and the transmission curves in each broadband for actual observatories such as GALEX, SDSS or
    Herschel. Refer to the description of the Band class for more information.

    The table below lists the broad bands available at the time of writing. The first column lists
    the complete band name; the second column indicates the corresponding pivot wavelength (a
    characteristic wavelength of the band at which the mean specific luminosity can be easily
    converted between wavelength and frequency representations).

    To specify a band, it suffices to enter just two segments of its name (case insensitive) that
    uniquely identify the band, seperated by an underscore or a space. For example, to select the
    HERSCHEL_PACS_100 band, one could enter "Herschel 100", "PACS 100", or "HERSCHEL_PACS_100".

    | %Band name | Pivot wavelength (micron)
    |------------|--------------------------
    | 2MASS_2MASS_J      |   1.24
    | 2MASS_2MASS_H      |   1.65
    | 2MASS_2MASS_KS     |   2.16
    | ALMA_ALMA_10       |    350
    | ALMA_ALMA_9        |    456
    | ALMA_ALMA_8        |    690
    | ALMA_ALMA_7        |    938
    | ALMA_ALMA_6        |   1244
    | GALEX_GALEX_FUV    |   0.15
    | GALEX_GALEX_NUV    |   0.23
    | GENERIC_JOHNSON_U  |   0.35
    | GENERIC_JOHNSON_B  |   0.44
    | GENERIC_JOHNSON_V  |   0.55
    | GENERIC_JOHNSON_R  |   0.69
    | GENERIC_JOHNSON_I  |   0.87
    | GENERIC_JOHNSON_J  |   1.24
    | GENERIC_JOHNSON_M  |   5.01
    | HERSCHEL_PACS_70   |   70.8
    | HERSCHEL_PACS_100  |    101
    | HERSCHEL_PACS_160  |    162
    | HERSCHEL_SPIRE_250 |    253
    | HERSCHEL_SPIRE_350 |    354
    | HERSCHEL_SPIRE_500 |    515
    | IRAS_IRAS_12       |   11.4
    | IRAS_IRAS_25       |   23.6
    | IRAS_IRAS_60       |   60.4
    | IRAS_IRAS_100      |    101
    | JCMT_SCUBA2_450    |    449
    | JCMT_SCUBA2_850    |    854
    | PLANCK_HFI_857     |    352
    | PLANCK_HFI_545     |    545
    | PLANCK_HFI_353     |    839
    | PLANCK_HFI_217     |   1368
    | SLOAN_SDSS_U       |   0.36
    | SLOAN_SDSS_G       |   0.47
    | SLOAN_SDSS_R       |   0.62
    | SLOAN_SDSS_I       |   0.75
    | SLOAN_SDSS_Z       |   0.89
    | SPITZER_IRAC_I1    |   3.55
    | SPITZER_IRAC_I2    |   4.50
    | SPITZER_IRAC_I3    |   5.72
    | SPITZER_IRAC_I4    |   7.88
    | SPITZER_MIPS_24    |   23.8
    | SPITZER_MIPS_70    |     72
    | SPITZER_MIPS_160   |    156
    | SWIFT_UVOT_UVW2    |   0.21
    | SWIFT_UVOT_UVM2    |   0.22
    | SWIFT_UVOT_UVW1    |   0.26
    | SWIFT_UVOT_U       |   0.35
    | SWIFT_UVOT_B       |   0.43
    | SWIFT_UVOT_V       |   0.54
    | TNG_OIG_U          |   0.37
    | TNG_OIG_B          |   0.44
    | TNG_OIG_V          |   0.54
    | TNG_OIG_R          |   0.64
    | TNG_NICS_J         |   1.28
    | TNG_NICS_H         |   1.63
    | TNG_NICS_K         |    2.2
    | UKIRT_UKIDSS_Z     |   0.88
    | UKIRT_UKIDSS_Y     |   1.03
    | UKIRT_UKIDSS_J     |   1.25
    | UKIRT_UKIDSS_H     |   1.64
    | UKIRT_UKIDSS_K     |   2.21
    | WISE_WISE_W1       |   3.39
    | WISE_WISE_W2       |   4.64
    | WISE_WISE_W3       |   12.6
    | WISE_WISE_W4       |   22.3

*/
class BroadBand : public Band
{
    ITEM_CONCRETE(BroadBand, Band, "a standard built-in broadband (transmission curve)")

        PROPERTY_STRING(bandName, "the name of the standard broadband (e.g. 'Johnson V' or 'PACS 100')")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked by classes that wish to hard-code the creation of a new
        BroadBand object (as opposed to creation through the ski file). Before the constructor
        returns, the newly created object is hooked up as a child to the specified parent in the
        simulation hierarchy (so it will automatically be deleted), and its setup() function has
        been called. The second argument specifies the name of a standard built-in broadband in the
        same format as it would be entered by a user in the configurable bandName property. */
    explicit BroadBand(SimulationItem* parent, string bandName);

protected:
    /** This function locates and opens the resource file corresponding to the configured band
        name. */
    void setupSelfBefore() override;

    //============= Functions required by base class =============

protected:
    /** This function returns the number of elements in the wavelength and transmission data arrays
        held by this subclass. */
    size_t dataSize() const override;

    /** This function returns a pointer to the first wavelength in the corresponding array held by
        this subclass. The number of elements in this array can be obtained through the dataSize()
        function. */
    const double* wavelengthData() const override;

    /** This function returns a pointer to the first transmission value in the corresponding array
        held by this subclass. The number of elements in this array can be obtained through the
        dataSize() function. The transmission values are normalized as described in the Band class
        header. */
    const double* transmissionData() const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    StoredTable<1> _table;
};

////////////////////////////////////////////////////////////////////

#endif
