/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TRUSTGRAPHITEGRAINCOMPOSITION_HPP
#define TRUSTGRAPHITEGRAINCOMPOSITION_HPP

#include "GrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The TrustGraphiteGrainComposition class represents the optical and calorimetric properties of
    graphite dust grains according to the dust model used for the TRUST benchmark (Gordon et al.
    2017, A&A, 603, 114) and the SHG benchmark (Camps et al. 2015, A&A 580, A87).

    The underlying data is provided by Karel Misselt as part of a download from the TRUST web site
    (http://ipag.osug.fr/RT13/RTTRUST/opa.php) and generally describes the BARE-GR-S model of
    Zubko, Dwek, and Arendt 2004, ApJS, 152, 211. The bulk mass density for the graphite grains is
    set to the value of 2240 kg/m3. */
class TrustGraphiteGrainComposition : public GrainComposition
{
    ITEM_CONCRETE(TrustGraphiteGrainComposition, GrainComposition, "a TRUST benchmark graphite dust grain composition")
    ITEM_END()

public:
    /** This constructor can be invoked by classes that wish to hard-code the creation of a new
        grain composition object of this type (as opposed to creation through the ski file). Before
        the constructor returns, the newly created object is hooked up as a child to the specified
        parent in the simulation hierarchy (so it will automatically be deleted), and it's setup()
        function has been called. */
    explicit TrustGraphiteGrainComposition(SimulationItem* parent);

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
