/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BEGEMANNPOROUSALUMINAGRAINCOMPOSITION_HPP
#define BEGEMANNPOROUSALUMINAGRAINCOMPOSITION_HPP

#include "GrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The BegemannPorousAluminaGrainComposition class represents the optical and calorimetric
    properties of porous alumina dust grains. The optical properties are taken from Begemann et al.
    (1997) for wavelengths >7.8 microns. The optical properties for wavelengths <7.8 microns are
    taken from Laor, A, & Draine, B.T. 1993 ApJ 402,441. The calorimetric properties are those for
    astronomical silicates given by Draine & Li (2007), originally from Draine & Lee (1984). The
    bulk mass density is set to 4020 kg/m3. */
class BegemannPorousAluminaGrainComposition : public GrainComposition
{
    ITEM_CONCRETE(BegemannPorousAluminaGrainComposition, GrainComposition,
                  "a Begemann porous alumina dust grain composition")
    ITEM_END()

public:
    /** This constructor can be invoked by classes that wish to hard-code the creation of a new
        grain composition object of this type (as opposed to creation through the ski file). Before
        the constructor returns, the newly created object is hooked up as a child to the specified
        parent in the simulation hierarchy (so it will automatically be deleted), and its setup()
        function has been called. */
    explicit BegemannPorousAluminaGrainComposition(SimulationItem* parent);

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
