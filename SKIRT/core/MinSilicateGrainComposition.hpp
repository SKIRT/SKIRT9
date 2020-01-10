/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MINSILICATEGRAINCOMPOSITION_HPP
#define MINSILICATEGRAINCOMPOSITION_HPP

#include "GrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The MinSilicateGrainComposition class represents the optical properties of amorphous silicate
    dust grains taken from Min et al. (2007, A&A, 462, 667). This model was designed to match the
    observations of interstellar dust towards the galactic center. The calorimetric properties are
    taken from the data for astronomical silicates supplied with the DustEM code, described by
    Compiègne et al. 2011 (AA, 525, A103), downloaded from http://www.ias.u-psud.fr/DUSTEM/. The
    bulk mass density is set to the value of 3090 kg/m3 specified by Min for silicate grains. */
class MinSilicateGrainComposition : public GrainComposition
{
    ITEM_CONCRETE(MinSilicateGrainComposition, GrainComposition, "a Min 2007 amorphous silicate dust grain composition")
    ITEM_END()

public:
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
