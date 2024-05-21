/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DRAINELI2007GRAPHITEGRAINCOMPOSITION_HPP
#define DRAINELI2007GRAPHITEGRAINCOMPOSITION_HPP

#include "GrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The DraineLi2007GraphiteGrainComposition class represents the optical and calorimetric properties of
    graphite dust grains. The optical properties are taken from Bruce Draine's website
    (https://www.astro.princeton.edu/~draine/dust/dust.diel.html) as those for graphite (Draine &
    Lee 1984; Laor & Draine 1993). The calorimetric properties follow the analytical enthalpy
    prescription for graphite given by equation (11) of Draine & Li (2007), properly integrated to
    obtain internal energy rather than heat capacity. The bulk mass density is set to the standard
    value of 2240 kg/m3 for graphite grains. */
class DraineLi2007GraphiteGrainComposition : public GrainComposition
{
    ITEM_CONCRETE(DraineLi2007GraphiteGrainComposition, GrainComposition, "a Draine graphite dust grain composition")
    ITEM_END()

public:
    /** This constructor can be invoked by classes that wish to hard-code the creation of a new
        grain composition object of this type (as opposed to creation through the ski file). Before
        the constructor returns, the newly created object is hooked up as a child to the specified
        parent in the simulation hierarchy (so it will automatically be deleted), and its setup()
        function has been called. */
    explicit DraineLi2007GraphiteGrainComposition(SimulationItem* parent);

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
