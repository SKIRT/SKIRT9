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
    | 2MASS_2MASS_J      | 1.2393
    | 2MASS_2MASS_H      | 1.6494
    | 2MASS_2MASS_KS     | 2.1638
    | ALMA_ALMA_10       | 349.89
    | ALMA_ALMA_9        | 456.2
    | ALMA_ALMA_8        | 689.59
    | ALMA_ALMA_7        | 937.98
    | ALMA_ALMA_6        | 1244.4
    | ALMA_ALMA_5        | 1616
    | ALMA_ALMA_4        | 2100.2
    | ALMA_ALMA_3        | 3043.4
    | EUCLID_VIS_VIS     | 0.71032
    | EUCLID_NISP_Y      | 1.0808
    | EUCLID_NISP_J      | 1.3644
    | EUCLID_NISP_H      | 1.7696
    | GALEX_GALEX_FUV    | 0.15351
    | GALEX_GALEX_NUV    | 0.23008
    | GENERIC_JOHNSON_U  | 0.35236
    | GENERIC_JOHNSON_B  | 0.44146
    | GENERIC_JOHNSON_V  | 0.55223
    | GENERIC_JOHNSON_R  | 0.68967
    | GENERIC_JOHNSON_I  | 0.87374
    | GENERIC_JOHNSON_J  | 1.2429
    | GENERIC_JOHNSON_M  | 5.0114
    | HERSCHEL_PACS_70   | 70.77
    | HERSCHEL_PACS_100  | 100.8
    | HERSCHEL_PACS_160  | 161.89
    | HERSCHEL_SPIRE_250 | 252.55
    | HERSCHEL_SPIRE_350 | 354.27
    | HERSCHEL_SPIRE_500 | 515.36
    | IRAS_IRAS_12       | 11.4
    | IRAS_IRAS_25       | 23.605
    | IRAS_IRAS_60       | 60.344
    | IRAS_IRAS_100      | 101.05
    | JCMT_SCUBA2_450    | 449.3
    | JCMT_SCUBA2_850    | 853.81
    | PLANCK_HFI_857     | 352.42
    | PLANCK_HFI_545     | 545.55
    | PLANCK_HFI_353     | 839.3
    | PLANCK_HFI_217     | 1367.6
    | PLANCK_HFI_143     | 2130.7
    | PLANCK_HFI_100     | 3001.1
    | PLANCK_LFI_70      | 4303
    | PLANCK_LFI_44      | 6845.9
    | PLANCK_LFI_30      | 10674
    | RUBIN_LSST_U       | 0.368
    | RUBIN_LSST_G       | 0.47823
    | RUBIN_LSST_R       | 0.62178
    | RUBIN_LSST_I       | 0.75323
    | RUBIN_LSST_Z       | 0.86851
    | RUBIN_LSST_Y       | 0.97301
    | SLOAN_SDSS_U       | 0.35565
    | SLOAN_SDSS_G       | 0.47024
    | SLOAN_SDSS_R       | 0.61755
    | SLOAN_SDSS_I       | 0.74899
    | SLOAN_SDSS_Z       | 0.89467
    | SPITZER_IRAC_I1    | 3.5508
    | SPITZER_IRAC_I2    | 4.496
    | SPITZER_IRAC_I3    | 5.7245
    | SPITZER_IRAC_I4    | 7.8842
    | SPITZER_MIPS_24    | 23.759
    | SPITZER_MIPS_70    | 71.987
    | SPITZER_MIPS_160   | 156.43
    | SWIFT_UVOT_UVW2    | 0.20551
    | SWIFT_UVOT_UVM2    | 0.22462
    | SWIFT_UVOT_UVW1    | 0.25804
    | SWIFT_UVOT_U       | 0.34628
    | SWIFT_UVOT_B       | 0.43496
    | SWIFT_UVOT_V       | 0.54254
    | TNG_OIG_U          | 0.37335
    | TNG_OIG_B          | 0.43975
    | TNG_OIG_V          | 0.53727
    | TNG_OIG_R          | 0.63917
    | TNG_NICS_J         | 1.2758
    | TNG_NICS_H         | 1.6265
    | TNG_NICS_K         | 2.2016
    | UKIRT_UKIDSS_Z     | 0.88263
    | UKIRT_UKIDSS_Y     | 1.0314
    | UKIRT_UKIDSS_J     | 1.2501
    | UKIRT_UKIDSS_H     | 1.6354
    | UKIRT_UKIDSS_K     | 2.2058
    | WISE_WISE_W1       | 3.3897
    | WISE_WISE_W2       | 4.6406
    | WISE_WISE_W3       | 12.568
    | WISE_WISE_W4       | 22.314

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
