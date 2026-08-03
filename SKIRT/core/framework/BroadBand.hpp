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

    The table below lists the broad bands included with SKIRT by default. The first column lists
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


    <b>James Webb Space Telescope</b>

    These broadband definitions are available only if the optional resource pack "ExtraBands" is
    installed (run './downloadResources.sh', see the installation guide for more info).

    | %Band name | Pivot wavelength (micron)
    |--------------------|------------------
    | JWST_NIRCAM_F070W  | 0.7039
    | JWST_NIRCAM_F090W  | 0.9022
    | JWST_NIRCAM_F115W  | 1.154
    | JWST_NIRCAM_F140M  | 1.405
    | JWST_NIRCAM_F150W  | 1.501
    | JWST_NIRCAM_F162M  | 1.627
    | JWST_NIRCAM_F164N  | 1.645
    | JWST_NIRCAM_F150W2 | 1.659
    | JWST_NIRCAM_F182M  | 1.845
    | JWST_NIRCAM_F187N  | 1.874
    | JWST_NIRCAM_F200W  | 1.989
    | JWST_NIRCAM_F210M  | 2.095
    | JWST_NIRCAM_F212N  | 2.121
    | JWST_NIRCAM_F250M  | 2.503
    | JWST_NIRCAM_F277W  | 2.762
    | JWST_NIRCAM_F300M  | 2.989
    | JWST_NIRCAM_F322W2 | 3.232
    | JWST_NIRCAM_F323N  | 3.237
    | JWST_NIRCAM_F335M  | 3.362
    | JWST_NIRCAM_F356W  | 3.568
    | JWST_NIRCAM_F360M  | 3.624
    | JWST_NIRCAM_F405N  | 4.052
    | JWST_NIRCAM_F410M  | 4.082
    | JWST_NIRCAM_F430M  | 4.281
    | JWST_NIRCAM_F444W  | 4.404
    | JWST_NIRCAM_F460M  | 4.630
    | JWST_NIRCAM_F466N  | 4.654
    | JWST_NIRCAM_F470N  | 4.708
    | JWST_NIRCAM_F480M  | 4.818
    | JWST_MIRI_F560W    | 5.635
    | JWST_MIRI_F770W    | 7.639
    | JWST_MIRI_F1000W   | 9.953
    | JWST_MIRI_F1130W   | 11.31
    | JWST_MIRI_F1280W   | 12.81
    | JWST_MIRI_F1500W   | 15.06
    | JWST_MIRI_F1800W   | 17.98
    | JWST_MIRI_F2100W   | 20.80
    | JWST_MIRI_F2550W   | 25.36

    <b>Hubble Space Telescope</b>

    These broadband definitions are available only if the optional resource pack "ExtraBands" is
    installed (run './downloadResources.sh', see the installation guide for more info).

    | %Band name | Pivot wavelength (micron)
    |--------------------------------|--------
    | HST_ACS_HRC_F220W              | 0.2254
    | HST_ACS_HRC_F250W              | 0.2716
    | HST_ACS_HRC_F330W              | 0.3363
    | HST_ACS_HRC_F344N              | 0.3434
    | HST_ACS_HRC_F435W              | 0.4323
    | HST_ACS_HRC_F475W              | 0.4776
    | HST_ACS_HRC_F502N              | 0.5023
    | HST_ACS_HRC_F550M              | 0.5580
    | HST_ACS_HRC_F555W              | 0.5356
    | HST_ACS_HRC_F606W              | 0.5887
    | HST_ACS_HRC_F625W              | 0.6295
    | HST_ACS_HRC_F658N              | 0.6584
    | HST_ACS_HRC_F660N              | 0.6599
    | HST_ACS_HRC_F775W              | 0.7665
    | HST_ACS_HRC_F814W              | 0.8100
    | HST_ACS_HRC_F850LP             | 0.9144
    | HST_ACS_HRC_F892N              | 0.8916
    | HST_ACS_HRC_POL_UV             | 0.6232
    | HST_ACS_HRC_POL_V              | 0.6959
    | HST_ACS_WFC_F435W              | 0.4330
    | HST_ACS_WFC_F475W              | 0.4747
    | HST_ACS_WFC_F502N              | 0.5023
    | HST_ACS_WFC_F550M              | 0.5581
    | HST_ACS_WFC_F555W              | 0.5361
    | HST_ACS_WFC_F606W              | 0.5922
    | HST_ACS_WFC_F625W              | 0.6312
    | HST_ACS_WFC_F658N              | 0.6584
    | HST_ACS_WFC_F660N              | 0.6599
    | HST_ACS_WFC_F775W              | 0.7693
    | HST_ACS_WFC_F814W              | 0.8046
    | HST_ACS_WFC_F850LP             | 0.9031
    | HST_ACS_WFC_F892N              | 0.8915
    | HST_ACS_WFC_POL_UV             | 0.6634
    | HST_ACS_WFC_POL_V              | 0.6911
    | HST_NICMOS1_F090M              | 0.9035
    | HST_NICMOS1_F095N              | 0.9537
    | HST_NICMOS1_F097N              | 0.9717
    | HST_NICMOS1_F108N              | 1.082
    | HST_NICMOS1_F110M              | 1.102
    | HST_NICMOS1_F110W              | 1.123
    | HST_NICMOS1_F113N              | 1.130
    | HST_NICMOS1_F140W              | 1.427
    | HST_NICMOS1_F145M              | 1.455
    | HST_NICMOS1_F160W              | 1.604
    | HST_NICMOS1_F164N              | 1.646
    | HST_NICMOS1_F165M              | 1.648
    | HST_NICMOS1_F166N              | 1.661
    | HST_NICMOS1_F170M              | 1.706
    | HST_NICMOS1_F187N              | 1.875
    | HST_NICMOS1_F190N              | 1.899
    | HST_NICMOS2_F110W              | 1.123
    | HST_NICMOS2_F160W              | 1.603
    | HST_NICMOS2_F165M              | 1.651
    | HST_NICMOS2_F171M              | 1.721
    | HST_NICMOS2_F180M              | 1.797
    | HST_NICMOS2_F187N              | 1.874
    | HST_NICMOS2_F187W              | 1.871
    | HST_NICMOS2_F190N              | 1.900
    | HST_NICMOS2_F204M              | 2.035
    | HST_NICMOS2_F205W              | 2.064
    | HST_NICMOS2_F207M              | 2.082
    | HST_NICMOS2_F212N              | 2.122
    | HST_NICMOS2_F215N              | 2.148
    | HST_NICMOS2_F216N              | 2.164
    | HST_NICMOS2_F222M              | 2.218
    | HST_NICMOS2_F237M              | 2.369
    | HST_NICMOS3_F108N              | 1.080
    | HST_NICMOS3_F110W              | 1.120
    | HST_NICMOS3_F113N              | 1.128
    | HST_NICMOS3_F150W              | 1.535
    | HST_NICMOS3_F160W              | 1.604
    | HST_NICMOS3_F164N              | 1.646
    | HST_NICMOS3_F166N              | 1.658
    | HST_NICMOS3_F175W              | 1.809
    | HST_NICMOS3_F187N              | 1.875
    | HST_NICMOS3_F190N              | 1.900
    | HST_NICMOS3_F196N              | 1.964
    | HST_NICMOS3_F200N              | 1.998
    | HST_NICMOS3_F212N              | 2.122
    | HST_NICMOS3_F215N              | 2.148
    | HST_NICMOS3_F222M              | 2.218
    | HST_NICMOS3_F240M              | 2.396
    | HST_NICMOS3_G096               | 1.019
    | HST_NICMOS3_G141               | 1.539
    | HST_NICMOS3_G206               | 2.084
    | HST_STIS_CCD_50CCD             | 0.5743
    | HST_STIS_CCD_50CCD_G230LB      | 0.2583
    | HST_STIS_CCD_50CCD_G230MB      | 0.2639
    | HST_STIS_CCD_50CCD_G430L       | 0.4516
    | HST_STIS_CCD_50CCD_G430M       | 0.4442
    | HST_STIS_CCD_50CCD_G750L       | 0.7151
    | HST_STIS_CCD_50CCD_G750M       | 0.7278
    | HST_STIS_CCD_50CORON           | 0.5743
    | HST_STIS_CCD_50CORON_G230LB    | 0.2583
    | HST_STIS_CCD_50CORON_G230MB    | 0.2639
    | HST_STIS_CCD_50CORON_G430L     | 0.4516
    | HST_STIS_CCD_50CORON_G430M     | 0.4442
    | HST_STIS_CCD_50CORON_G750L     | 0.7151
    | HST_STIS_CCD_50CORON_G750M     | 0.7278
    | HST_STIS_CCD_F28X50LP          | 0.7200
    | HST_STIS_CCD_F28X50LP_G230LB   | 0.2583
    | HST_STIS_CCD_F28X50LP_G230MB   | 0.2639
    | HST_STIS_CCD_F28X50LP_G430L    | 0.5600
    | HST_STIS_CCD_F28X50LP_G430M    | 0.5527
    | HST_STIS_CCD_F28X50LP_G750L    | 0.7283
    | HST_STIS_CCD_F28X50LP_G750M    | 0.7348
    | HST_STIS_CCD_F28X50OII         | 0.3738
    | HST_STIS_CCD_F28X50OIII        | 0.5006
    | HST_STIS_CCD_F28X50OIII_G430L  | 0.5006
    | HST_STIS_CCD_F28X50OIII_G430M  | 0.5006
    | HST_STIS_CCD_F28X50OII_G430L   | 0.3738
    | HST_STIS_CCD_F28X50OII_G430M   | 0.3738
    | HST_STIS_FUV_25MAMA            | 0.1369
    | HST_STIS_FUV_25MAMA_G140L      | 0.1353
    | HST_STIS_FUV_25MAMA_G140M      | 0.1349
    | HST_STIS_FUV_F25LYA            | 0.1243
    | HST_STIS_FUV_F25LYA_G140L      | 0.1246
    | HST_STIS_FUV_F25LYA_G140M      | 0.1246
    | HST_STIS_FUV_F25ND3            | 0.1371
    | HST_STIS_FUV_F25ND3_G140L      | 0.1356
    | HST_STIS_FUV_F25ND3_G140M      | 0.1352
    | HST_STIS_FUV_F25ND5            | 0.1379
    | HST_STIS_FUV_F25ND5_G140L      | 0.1361
    | HST_STIS_FUV_F25ND5_G140M      | 0.1358
    | HST_STIS_FUV_F25NDQ1           | 0.1409
    | HST_STIS_FUV_F25NDQ1_G140L     | 0.1386
    | HST_STIS_FUV_F25NDQ1_G140M     | 0.1382
    | HST_STIS_FUV_F25NDQ2           | 0.1398
    | HST_STIS_FUV_F25NDQ2_G140L     | 0.1376
    | HST_STIS_FUV_F25NDQ2_G140M     | 0.1372
    | HST_STIS_FUV_F25NDQ3           | 0.1401
    | HST_STIS_FUV_F25NDQ3_G140L     | 0.1380
    | HST_STIS_FUV_F25NDQ3_G140M     | 0.1376
    | HST_STIS_FUV_F25NDQ4           | 0.1427
    | HST_STIS_FUV_F25NDQ4_G140L     | 0.1401
    | HST_STIS_FUV_F25NDQ4_G140M     | 0.1397
    | HST_STIS_FUV_F25QTZ            | 0.1595
    | HST_STIS_FUV_F25QTZ_G140L      | 0.1561
    | HST_STIS_FUV_F25QTZ_G140M      | 0.1564
    | HST_STIS_FUV_F25SRF2           | 0.1453
    | HST_STIS_FUV_F25SRF2_G140L     | 0.1426
    | HST_STIS_FUV_F25SRF2_G140M     | 0.1423
    | HST_STIS_NUV_25MAMA            | 0.2260
    | HST_STIS_NUV_25MAMA_G230L      | 0.2418
    | HST_STIS_NUV_25MAMA_G230M      | 0.2431
    | HST_STIS_NUV_25MAMA_PRISM      | 0.2280
    | HST_STIS_NUV_F25CIII           | 0.2010
    | HST_STIS_NUV_F25CIII_G230L     | 0.2008
    | HST_STIS_NUV_F25CIII_G230M     | 0.2012
    | HST_STIS_NUV_F25CIII_PRISM     | 0.2015
    | HST_STIS_NUV_F25CN182          | 0.2003
    | HST_STIS_NUV_F25CN182_G230L    | 0.2081
    | HST_STIS_NUV_F25CN182_G230M    | 0.2102
    | HST_STIS_NUV_F25CN182_PRISM    | 0.2014
    | HST_STIS_NUV_F25CN270          | 0.2711
    | HST_STIS_NUV_F25CN270_G230L    | 0.2710
    | HST_STIS_NUV_F25CN270_G230M    | 0.2715
    | HST_STIS_NUV_F25CN270_PRISM    | 0.2712
    | HST_STIS_NUV_F25MGII           | 0.2802
    | HST_STIS_NUV_F25MGII_G230L     | 0.2802
    | HST_STIS_NUV_F25MGII_G230M     | 0.2802
    | HST_STIS_NUV_F25MGII_PRISM     | 0.2802
    | HST_STIS_NUV_F25ND3            | 0.2356
    | HST_STIS_NUV_F25ND3_G230L      | 0.2517
    | HST_STIS_NUV_F25ND3_G230M      | 0.2528
    | HST_STIS_NUV_F25ND3_PRISM      | 0.2382
    | HST_STIS_NUV_F25ND5            | 0.2630
    | HST_STIS_NUV_F25ND5_G230L      | 0.2730
    | HST_STIS_NUV_F25ND5_G230M      | 0.2732
    | HST_STIS_NUV_F25ND5_PRISM      | 0.2663
    | HST_STIS_NUV_F25NDQ1           | 0.2305
    | HST_STIS_NUV_F25NDQ1_G230L     | 0.2441
    | HST_STIS_NUV_F25NDQ1_G230M     | 0.2454
    | HST_STIS_NUV_F25NDQ1_PRISM     | 0.2336
    | HST_STIS_NUV_F25NDQ2           | 0.2335
    | HST_STIS_NUV_F25NDQ2_G230L     | 0.2480
    | HST_STIS_NUV_F25NDQ2_G230M     | 0.2494
    | HST_STIS_NUV_F25NDQ2_PRISM     | 0.2366
    | HST_STIS_NUV_F25NDQ3           | 0.2432
    | HST_STIS_NUV_F25NDQ3_G230L     | 0.2550
    | HST_STIS_NUV_F25NDQ3_G230M     | 0.2561
    | HST_STIS_NUV_F25NDQ3_PRISM     | 0.2465
    | HST_STIS_NUV_F25NDQ4           | 0.2500
    | HST_STIS_NUV_F25NDQ4_G230L     | 0.2605
    | HST_STIS_NUV_F25NDQ4_G230M     | 0.2613
    | HST_STIS_NUV_F25NDQ4_PRISM     | 0.2542
    | HST_STIS_NUV_F25QTZ            | 0.2358
    | HST_STIS_NUV_F25QTZ_G230L      | 0.2425
    | HST_STIS_NUV_F25QTZ_G230M      | 0.2437
    | HST_STIS_NUV_F25QTZ_PRISM      | 0.2390
    | HST_STIS_NUV_F25SRF2           | 0.2302
    | HST_STIS_NUV_F25SRF2_G230L     | 0.2426
    | HST_STIS_NUV_F25SRF2_G230M     | 0.2438
    | HST_STIS_NUV_F25SRF2_PRISM     | 0.2339
    | HST_WFC3_IR_F098M              | 0.9863
    | HST_WFC3_IR_F105W              | 1.055
    | HST_WFC3_IR_F110W              | 1.153
    | HST_WFC3_IR_F125W              | 1.249
    | HST_WFC3_IR_F126N              | 1.259
    | HST_WFC3_IR_F127M              | 1.274
    | HST_WFC3_IR_F128N              | 1.284
    | HST_WFC3_IR_F130N              | 1.301
    | HST_WFC3_IR_F132N              | 1.319
    | HST_WFC3_IR_F139M              | 1.384
    | HST_WFC3_IR_F140W              | 1.392
    | HST_WFC3_IR_F153M              | 1.533
    | HST_WFC3_IR_F160W              | 1.537
    | HST_WFC3_IR_F164N              | 1.645
    | HST_WFC3_IR_F167N              | 1.668
    | HST_WFC3_IR_G102               | 0.9990
    | HST_WFC3_IR_G141               | 1.389
    | HST_WFC3_UVIS1_F200LP          | 0.4972
    | HST_WFC3_UVIS1_F218W           | 0.2225
    | HST_WFC3_UVIS1_F225W           | 0.2371
    | HST_WFC3_UVIS1_F275W           | 0.2709
    | HST_WFC3_UVIS1_F280N           | 0.2797
    | HST_WFC3_UVIS1_F300X           | 0.2820
    | HST_WFC3_UVIS1_F336W           | 0.3354
    | HST_WFC3_UVIS1_F343N           | 0.3435
    | HST_WFC3_UVIS1_F350LP          | 0.5874
    | HST_WFC3_UVIS1_F373N           | 0.3730
    | HST_WFC3_UVIS1_F390M           | 0.3897
    | HST_WFC3_UVIS1_F390W           | 0.3924
    | HST_WFC3_UVIS1_F395N           | 0.3955
    | HST_WFC3_UVIS1_F410M           | 0.4109
    | HST_WFC3_UVIS1_F438W           | 0.4326
    | HST_WFC3_UVIS1_F467M           | 0.4682
    | HST_WFC3_UVIS1_F469N           | 0.4688
    | HST_WFC3_UVIS1_F475W           | 0.4773
    | HST_WFC3_UVIS1_F475X           | 0.4941
    | HST_WFC3_UVIS1_F487N           | 0.4871
    | HST_WFC3_UVIS1_F502N           | 0.5010
    | HST_WFC3_UVIS1_F547M           | 0.5448
    | HST_WFC3_UVIS1_F555W           | 0.5308
    | HST_WFC3_UVIS1_F600LP          | 0.7468
    | HST_WFC3_UVIS1_F606W           | 0.5889
    | HST_WFC3_UVIS1_F621M           | 0.6219
    | HST_WFC3_UVIS1_F625W           | 0.6243
    | HST_WFC3_UVIS1_F631N           | 0.6304
    | HST_WFC3_UVIS1_F645N           | 0.6453
    | HST_WFC3_UVIS1_F656N           | 0.6562
    | HST_WFC3_UVIS1_F657N           | 0.6567
    | HST_WFC3_UVIS1_F658N           | 0.6586
    | HST_WFC3_UVIS1_F665N           | 0.6656
    | HST_WFC3_UVIS1_F673N           | 0.6766
    | HST_WFC3_UVIS1_F680N           | 0.6878
    | HST_WFC3_UVIS1_F689M           | 0.6877
    | HST_WFC3_UVIS1_F763M           | 0.7614
    | HST_WFC3_UVIS1_F775W           | 0.7651
    | HST_WFC3_UVIS1_F814W           | 0.8039
    | HST_WFC3_UVIS1_F845M           | 0.8439
    | HST_WFC3_UVIS1_F850LP          | 0.9176
    | HST_WFC3_UVIS1_F953N           | 0.9531
    | HST_WFC3_UVIS1_FQ232N          | 0.2327
    | HST_WFC3_UVIS1_FQ243N          | 0.2421
    | HST_WFC3_UVIS1_FQ378N          | 0.3792
    | HST_WFC3_UVIS1_FQ387N          | 0.3874
    | HST_WFC3_UVIS1_FQ422M          | 0.4219
    | HST_WFC3_UVIS1_FQ436N          | 0.4367
    | HST_WFC3_UVIS1_FQ437N          | 0.4371
    | HST_WFC3_UVIS1_FQ492N          | 0.4933
    | HST_WFC3_UVIS1_FQ508N          | 0.5091
    | HST_WFC3_UVIS1_FQ575N          | 0.5757
    | HST_WFC3_UVIS1_FQ619N          | 0.6198
    | HST_WFC3_UVIS1_FQ634N          | 0.6349
    | HST_WFC3_UVIS1_FQ672N          | 0.6717
    | HST_WFC3_UVIS1_FQ674N          | 0.6731
    | HST_WFC3_UVIS1_FQ727N          | 0.7276
    | HST_WFC3_UVIS1_FQ750N          | 0.7502
    | HST_WFC3_UVIS1_FQ889N          | 0.8892
    | HST_WFC3_UVIS1_FQ906N          | 0.9058
    | HST_WFC3_UVIS1_FQ924N          | 0.9248
    | HST_WFC3_UVIS1_FQ937N          | 0.9373
    | HST_WFC3_UVIS1_G280            | 0.3749
    | HST_WFC3_UVIS2_F200LP          | 0.4875
    | HST_WFC3_UVIS2_F218W           | 0.2221
    | HST_WFC3_UVIS2_F225W           | 0.2358
    | HST_WFC3_UVIS2_F275W           | 0.2703
    | HST_WFC3_UVIS2_F280N           | 0.2797
    | HST_WFC3_UVIS2_F300X           | 0.2805
    | HST_WFC3_UVIS2_F336W           | 0.3355
    | HST_WFC3_UVIS2_F343N           | 0.3435
    | HST_WFC3_UVIS2_F350LP          | 0.5851
    | HST_WFC3_UVIS2_F373N           | 0.3730
    | HST_WFC3_UVIS2_F390M           | 0.3897
    | HST_WFC3_UVIS2_F390W           | 0.3921
    | HST_WFC3_UVIS2_F395N           | 0.3955
    | HST_WFC3_UVIS2_F410M           | 0.4109
    | HST_WFC3_UVIS2_F438W           | 0.4325
    | HST_WFC3_UVIS2_F467M           | 0.4682
    | HST_WFC3_UVIS2_F469N           | 0.4688
    | HST_WFC3_UVIS2_F475W           | 0.4772
    | HST_WFC3_UVIS2_F475X           | 0.4937
    | HST_WFC3_UVIS2_F487N           | 0.4871
    | HST_WFC3_UVIS2_F502N           | 0.5010
    | HST_WFC3_UVIS2_F547M           | 0.5447
    | HST_WFC3_UVIS2_F555W           | 0.5308
    | HST_WFC3_UVIS2_F600LP          | 0.7454
    | HST_WFC3_UVIS2_F606W           | 0.5888
    | HST_WFC3_UVIS2_F621M           | 0.6219
    | HST_WFC3_UVIS2_F625W           | 0.6242
    | HST_WFC3_UVIS2_F631N           | 0.6304
    | HST_WFC3_UVIS2_F645N           | 0.6453
    | HST_WFC3_UVIS2_F656N           | 0.6562
    | HST_WFC3_UVIS2_F657N           | 0.6567
    | HST_WFC3_UVIS2_F658N           | 0.6586
    | HST_WFC3_UVIS2_F665N           | 0.6656
    | HST_WFC3_UVIS2_F673N           | 0.6766
    | HST_WFC3_UVIS2_F680N           | 0.6877
    | HST_WFC3_UVIS2_F689M           | 0.6877
    | HST_WFC3_UVIS2_F763M           | 0.7613
    | HST_WFC3_UVIS2_F775W           | 0.7648
    | HST_WFC3_UVIS2_F814W           | 0.8029
    | HST_WFC3_UVIS2_F845M           | 0.8437
    | HST_WFC3_UVIS2_F850LP          | 0.9170
    | HST_WFC3_UVIS2_F953N           | 0.9531
    | HST_WFC3_UVIS2_FQ232N          | 0.2327
    | HST_WFC3_UVIS2_FQ243N          | 0.2421
    | HST_WFC3_UVIS2_FQ378N          | 0.3792
    | HST_WFC3_UVIS2_FQ387N          | 0.3874
    | HST_WFC3_UVIS2_FQ422M          | 0.4219
    | HST_WFC3_UVIS2_FQ436N          | 0.4367
    | HST_WFC3_UVIS2_FQ437N          | 0.4371
    | HST_WFC3_UVIS2_FQ492N          | 0.4933
    | HST_WFC3_UVIS2_FQ508N          | 0.5091
    | HST_WFC3_UVIS2_FQ575N          | 0.5757
    | HST_WFC3_UVIS2_FQ619N          | 0.6198
    | HST_WFC3_UVIS2_FQ634N          | 0.6349
    | HST_WFC3_UVIS2_FQ672N          | 0.6717
    | HST_WFC3_UVIS2_FQ674N          | 0.6731
    | HST_WFC3_UVIS2_FQ727N          | 0.7276
    | HST_WFC3_UVIS2_FQ750N          | 0.7502
    | HST_WFC3_UVIS2_FQ889N          | 0.8892
    | HST_WFC3_UVIS2_FQ906N          | 0.9058
    | HST_WFC3_UVIS2_FQ924N          | 0.9248
    | HST_WFC3_UVIS2_FQ937N          | 0.9373
    | HST_WFC3_UVIS2_G280            | 0.3652

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
