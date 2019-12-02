/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef POLARIZEDSILICATEGRAINCOMPOSITION_HPP
#define POLARIZEDSILICATEGRAINCOMPOSITION_HPP

#include "GrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The PolarizedSilicateGrainComposition class represents the optical and calorimetric properties
    of spherical silicate dust grains with support for polarization. The optical data, including
    scattering and absorption efficiency coefficients and Mueller matrix coefficients, are
    converted from a text file in the format as used by the STOKES code, a code for radiative
    transfer with polarization written by Rene Goosmann and Frédéric Marin. The data files were
    obtained from the STOKES team via Marko Stalevski <mstalevski@aob.rs>.

    The class also offers an option to specify the name of a stored table file providing the
    properties for polarized emission by spheroidal grains. This is a temporary provision; future
    implementations should provide a built-in resource.

    The calorimetric properties follow the analytical enthalpy prescription for silicate given by
    equation (9) of Draine & Li (2001), properly integrated to obtain internal energy rather than
    heat capacity. The bulk mass density is set to the standard value of 3000 kg/m3 for silicate
    grains. */
class PolarizedSilicateGrainComposition : public GrainComposition
{
    ITEM_CONCRETE(PolarizedSilicateGrainComposition, GrainComposition,
                  "a silicate dust grain composition with support for polarization")

    PROPERTY_STRING(spheroidalEmissionTable,
                    "the name of the file tabulating properties for polarized emission by spheroidal grains")
        ATTRIBUTE_REQUIRED_IF(spheroidalEmissionTable, "false")
        ATTRIBUTE_DISPLAYED_IF(spheroidalEmissionTable, "Level3")

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

    /** This function returns the name of the stored table resource tabulating the coefficients of
        the Mueller matrix as a function of wavelength, grain size, and scattering angle. */
    string resourceNameForMuellerMatrix() const override;

    /** This function returns the name of the stored table resource tabulating the properties
        relevant for polarized emission by spheroidal grains (absorption efficiencies and
        corresponding linear polarization efficiencies as a function of grain size, wavelength, and
        emission angle. The current implementation returns a string configured by the user; future
        implementations should return the name of a built-in resource. */
    virtual string resourceNameForSpheroidalEmission() const override;

    /** This function returns the name of the stored table resource tabulating the specific
        enthalpies per unit volume as a function of temperature. */
    string resourceNameForEnthalpies() const override;
};

////////////////////////////////////////////////////////////////////

#endif
