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

    The calorimetric properties follow the analytical enthalpy prescription for silicate given by
    equation (9) of Draine & Li (2001), properly integrated to obtain internal energy rather than
    heat capacity. The bulk mass density is set to the standard value of 3000 kg/m3 for silicate
    grains. */
class PolarizedSilicateGrainComposition : public GrainComposition
{
    ITEM_CONCRETE(PolarizedSilicateGrainComposition, GrainComposition,
                  "a silicate dust grain composition with support for polarization")
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

    /** This function returns the name of the stored table resource tabulating the specific
        enthalpies per unit volume as a function of temperature. */
    string resourceNameForEnthalpies() const override;
};

////////////////////////////////////////////////////////////////////

#endif
