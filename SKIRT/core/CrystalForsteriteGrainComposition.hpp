/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CRYSTALFORSTERITEGRAINCOMPOSITION_HPP
#define CRYSTALFORSTERITEGRAINCOMPOSITION_HPP

#include "GrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The CrystalForsteriteGrainComposition class represents the optical and calorimetric properties
    of crystalline silicate Mg2SiO4 grains, for which Michiel Min <michiel.min.science@gmail.com>
    prepared the data. The refractive index data was taken from Suto et al. 2006, using the lowest
    temperature 50K. Further data was obtained from the Jena group (Fabian 2001, Zeidler 2011) for
    UV to near-IR. Below 0.2 micron the results are extrapolated using theoretical formulas. The
    calculations were performed with DHS using \f$f_\mathrm{max}=0.8\f$ (see Min et al. 2005). The
    calorimetric properties are taken from the data for astronomical silicates supplied with the
    DustEM code, described by Compiègne et al. 2011 (AA, 525, A103), downloaded from
    http://www.ias.u-psud.fr/DUSTEM/. The bulk mass density is set to the value of 3330 kg/m3
    specified by Min. */
class CrystalForsteriteGrainComposition : public GrainComposition
{
    ITEM_CONCRETE(CrystalForsteriteGrainComposition, GrainComposition,
                  "a crystalline Forsterite dust grain composition")
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
