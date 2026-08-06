/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MIESILICATEGRAINCOMPOSITION_HPP
#define MIESILICATEGRAINCOMPOSITION_HPP

#include "GrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The MieSilicateGrainComposition class represents the optical and calorimetric properties of
    spherical silicate dust grains, based on the dielectric constants of astronomical silicate
    taken from Bruce Draine's website.

    The optical properties are obtained from a file containing precomputed values calculated using
    Mie theory. Sebastian Wolf's Mie program MieX was used to do the computations (it allows for
    large grains). This dust mixture is similar to the DraineSilicateGrainComposition class, except
    that the current class contains optical properties over a much wider grain size range: the
    current dust mixture contains 301 grain sizes ranging from 1 nm to 1 mm, whereas
    DraineSilicateGrainComposition contains 81 grain sizes from 1 nm to 10 micron.

    The calorimetric properties follow the analytical enthalpy prescription for silicate given by
    equation (11) of Draine & Li (2001), properly integrated to obtain internal energy rather than
    heat capacity. The bulk mass density is set to the standard value of 3000 kg/m3 for silicate
    grains. */
class MieSilicateGrainComposition : public GrainComposition
{
    ITEM_CONCRETE(MieSilicateGrainComposition, GrainComposition, "a MieX-based Draine silicate dust grain composition")
    ITEM_END()

public:
    /** This constructor can be invoked by  classes that wish to hard-code the creation of
        a new grain composition object of this type (as opposed to creation through the ski file).
        Before the constructor returns, the newly created object is hooked up as a child to the
        specified parent in the simulation hierarchy (so it will automatically be deleted), and
        it's setup() function has been called. */
    explicit MieSilicateGrainComposition(SimulationItem* parent);

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
