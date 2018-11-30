/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CRYSTALENSTATITEGRAINCOMPOSITION_HPP
#define CRYSTALENSTATITEGRAINCOMPOSITION_HPP

#include "GrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The CrystalEnstatiteGrainComposition class represents the optical and calorimetric properties
    of crystalline silicate MgSiO3 grains, for which Michiel Min <michiel.min.science@gmail.com>
    prepared the data. The refractive index data was taken from Jaeger et al. 1998. Further data
    was obtained from the Jena group (Fabian 2001, Zeidler 2011) for UV to near-IR. Below 0.2
    micron the results are extrapolated using theoretical formulas. The calculations were performed
    with DHS using \f$f_\mathrm{max}=0.8\f$ (see Min et al. 2005). The calorimetric properties are
    taken from the data for astronomical silicates supplied with the DustEM code, described by
    Compiègne et al. 2011 (AA, 525, A103), downloaded from http://www.ias.u-psud.fr/DUSTEM/. The
    bulk mass density is set to the value of 2800 kg/m3 specified by Min. */
class CrystalEnstatiteGrainComposition : public GrainComposition
{
    ITEM_CONCRETE(CrystalEnstatiteGrainComposition, GrainComposition, "a crystalline Enstatite dust grain composition")
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
