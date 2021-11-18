/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SINGLEWAVELENGTHSED_HPP
#define SINGLEWAVELENGTHSED_HPP

#include "LineSED.hpp"

////////////////////////////////////////////////////////////////////

/** The SingleWavelengthSED class implements a spectral energy distribution in the form of a
    Dirac-delta function, i.e. a single emission line with zero width. All photon packets are thus
    emitted at a single, configurable wavelength. The default value for this emission wavelength is
    the central wavelength of the hydrogen Lyman-alpha line.

    Because the specific luminosity of the spectrum is undefined, there are important restrictions
    on the use of this %SED. See the description of the LineSED class for more information. */
class SingleWavelengthSED : public LineSED
{
    ITEM_CONCRETE(SingleWavelengthSED, LineSED, "a single-wavelength SED in the form of a Dirac-delta function")

        PROPERTY_DOUBLE(wavelength, "the single emission wavelength")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")
        ATTRIBUTE_DEFAULT_VALUE(wavelength, "1215.67 Angstrom")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function return the single, configured emission wavelength and an arbitrary
        corresponding relative luminosity. */
    void getWavelengthsAndLuminosities(Array& lambdav, Array& Lv) const override;

    /** This function always returns the single emission wavelength. It is overridden here to avoid
        the random number generation needed in the base class for selecting between multiple lines.
        */
    double generateWavelength() const override;
};

////////////////////////////////////////////////////////////////////

#endif
