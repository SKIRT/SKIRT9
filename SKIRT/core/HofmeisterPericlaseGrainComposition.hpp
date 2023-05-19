/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef HOFMEISTERPERICLASEGRAINCOMPOSITION_HPP
#define HOFMEISTERPERICLASEGRAINCOMPOSITION_HPP

#include "GrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The HofmeisterPericlaseGrainComposition class represents the optical and calorimetric
    properties of periclase dust grains. The optical properties are taken from Hofmeister et al.
    (2003). The calorimetric properties are those for astronomical silicates given by Draine & Li
    (2007), originally from Draine & Lee (1984). The bulk mass density is set to 3560 kg/m3. The
    grain shape is assumed to be a Continuous Distribution of Ellipsoids (CDE) and the grain size
    corresponds to a volume that a spherical grain with the given radius would have. */
class HofmeisterPericlaseGrainComposition : public GrainComposition
{
    ITEM_CONCRETE(HofmeisterPericlaseGrainComposition, GrainComposition,
                  "a Hofmeister periclase dust grain composition")
    ITEM_END()

public:
    /** This constructor can be invoked by classes that wish to hard-code the creation of a new
        grain composition object of this type (as opposed to creation through the ski file). Before
        the constructor returns, the newly created object is hooked up as a child to the specified
        parent in the simulation hierarchy (so it will automatically be deleted), and its setup()
        function has been called. */
    explicit HofmeisterPericlaseGrainComposition(SimulationItem* parent);

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
