/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTEMGRAINCOMPOSITION_HPP
#define DUSTEMGRAINCOMPOSITION_HPP

#include "GrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The DustEmGrainComposition class represents the optical and calorimetric properties of one of a
    selection of dust grain types taken from the data supplied with the DustEM code, described by
    Compiègne et al. 2011 (AA, 525, A103), downloaded from http://www.ias.u-psud.fr/DUSTEM/.

    Through the \em grainType attribute, this class can be configured to represent any one of the
    following DustEM grain materials (text copied from the DustEM documentation). Both the optical
    and calorimetric properties are read from the DustEM data:

    - aSil: astronomical silicates. The Q-values have been generated with the Mie method. The
    refractive index and heat capacity are as in Draine & Li (2007), originally from Draine & Lee
    (1984).

    - Gra: graphite grains. The Q-values have been generated with the Mie method for spherical
    particles. The refractive index and heat capacity are as described in Draine & Li (2001) and Li
    & Draine (2001).

    - PAH0_DL07, PAH1_DL07: PAHs as defined in Draine & Li (2007) for neutral and singly charged
    PAHs, respectively. It takes into account the transition to graphite optical properties for a >
    5 nm.

    - PAH0_MC10, PAH1_MC10: PAHs as defined in Compiègne et al. (2011) for neutral and singly
    charged PAHs, respectively.

    - CM20: carbonaceous grains as defined in Jones et al. (2013). Grains up to 20 nm consist
    purely of aromatic-rich H-poor amorphous carbon, a-C, whereas bigger grains have a core/mantle
    structure, where the core is made of aliphatic-rich H-rich carbon, a-C:H, and the mantle of a-C
    with a thickness of 20 nm.

    - aOlM5 + aPyM5: silicates as defined in Koehler et al. (2014). These core/mantle grains
    consist of amorphous silicate cores (olivine and pyroxene-normative composition, Mg-rich) and 5
    nm mantles of a-C.

    The bulk mass density must be configured by the user, although a reasonable default value is
    provided for each of the supported grain types: 3000 kg/m3 for astronomical silicate, 2240
    kg/m3 for graphite and for the PAHs, 1600 kg/m3 for the H-rich carbonaceous grains, and 2190
    kg/m3 for the amorphous olivine and pyroxene.

    The last three listed grain types are used in the Themis dust model described by Jones et al.
    2017 (A&A, 602, A46) and the references therein. Specifically, CM20 represents the amorphous
    carbonaceous dust grains in that model, aOlM5 (olivine) represents the amorphous silicates with
    forsterite-normative composition, and aPyM5 (pyroxene) represents the amorphous silicates with
    enstatite-normative composition. */
class DustEmGrainComposition : public GrainComposition
{
    /** The enumeration type indicating the DustEM grain type represented by this class instance. */
    ENUM_DEF(GrainType, aSil, Gra, PAH0DL07, PAH1DL07, PAH0MC10, PAH1MC10, CM20, aOlM5, aPyM5)
        ENUM_VAL(GrainType, aSil, "astronomical silicate grains (Draine & Li 2007)")
        ENUM_VAL(GrainType, Gra, "graphite grains (Draine & Li 2001; Li & Draine 2001)")
        ENUM_VAL(GrainType, PAH0DL07, "neutral PAHs (Draine & Li 2007)")
        ENUM_VAL(GrainType, PAH1DL07, "ionized PAHs (Draine & Li 2007)")
        ENUM_VAL(GrainType, PAH0MC10, "neutral PAHs (Compiègne et al. 2011)")
        ENUM_VAL(GrainType, PAH1MC10, "ionized PAHs (Compiègne et al. 2011)")
        ENUM_VAL(GrainType, CM20, "amorphous hydro-carbon grains (Jones et al. 2013)")
        ENUM_VAL(GrainType, aOlM5, "amorphous olivine/forsterite grains (Koehler et al. 2014)")
        ENUM_VAL(GrainType, aPyM5, "amorphous pyroxene/enstatite grains (Koehler et al. 2014)")
    ENUM_END()

    ITEM_CONCRETE(DustEmGrainComposition, GrainComposition, "a dust grain composition based on DustEM data")

        PROPERTY_ENUM(grainType, GrainType, "the DustEM grain type")

        PROPERTY_DOUBLE(bulkMassDensity, "the bulk mass density for this grain material")
        ATTRIBUTE_QUANTITY(bulkMassDensity, "bulkmassdensity")
        ATTRIBUTE_MIN_VALUE(bulkMassDensity, "[100 kg/m3")
        ATTRIBUTE_MAX_VALUE(bulkMassDensity, "10000 kg/m3]")
        ATTRIBUTE_DEFAULT_VALUE(
            bulkMassDensity,
            "grainTypeaSil:3500 kg/m3;grainTypeCM20:1600 kg/m3;grainTypeaOlM5|grainTypeaPyM5:2190 kg/m3;2240 kg/m3")

    ITEM_END()

public:
    /** This constructor can be invoked by classes that wish to hard-code the creation of
        a new grain composition object of this type (as opposed to creation through the ski file).
        Before the constructor returns, the newly created object is hooked up as a child to the
        specified parent in the simulation hierarchy (so it will automatically be deleted), the
        grain type attribute has been set to the value specified in the second argument, the bulk
        mass density has been set to the value specified in the third argument, and the setup()
        function has been called. */
    explicit DustEmGrainComposition(SimulationItem* parent, GrainType grainType, double bulkMassDensity);

    /** This function returns a brief human-readable identifier for this grain composition. */
    string name() const override;

    /** This function returns the bulk mass density of this grain material. */
    double bulkDensity() const override;

    /** This function returns the name of the stored table resource tabulating the basic optical
        properties (absorption and scattering efficiencies and asymmetry parameter) as a function
        of wavelength and grain size. */
    string resourceNameForOpticalProps() const override;

    /** This function returns the name of the stored table resource tabulating the specific
        enthalpies per unit volume as a function of temperature. */
    string resourceNameForEnthalpies() const override;
};

////////////////////////////////////////////////////////////////////

#endif
