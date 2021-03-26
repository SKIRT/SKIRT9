/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LINESED_HPP
#define LINESED_HPP

#include "Array.hpp"
#include "SED.hpp"

////////////////////////////////////////////////////////////////////

/** LineSED is an abstract class for representing spectral energy distributions that consist of one
    or more discrete emission lines with zero width. Mathematically, the specific luminosity is
    infinite at the line wavelengths and zero everywhere else, while the integrated luminosity
    (over all lines in the source wavelength range) is still normalized to unity. As a result, the
    specific luminosity is numerically undefined. This imposes important restrictions on the use of the %SED:

    - A LineSED must be normalized by defining the integrated luminosity over a wavelength range
    that includes one or more of the line wavelengths. It \em cannot be normalized by defining the
    specific luminosity at a given wavelength or for a broadband.

    - The wavelength bias property of the corresponding source must be set to zero, because the
    mechanism for wavelength biasing relies on the specific luminosity.

    The subclass must provide the list of wavelengths and corresponding relative luminosities, and
    this abstract class handles everything else. */
class LineSED : public SED
{
    ITEM_ABSTRACT(LineSED, SED, "a spectral energy distribution consisting of discrete lines")
        ATTRIBUTE_TYPE_ALLOWED_IF(LineSED, "Panchromatic")
        ATTRIBUTE_TYPE_DISPLAYED_IF(LineSED, "Level2")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function asks the subclass to load the wavelength/luminosity pairs for one or more
        lines and precalculates the information needed by the other functions in this class. */
    void setupSelfBefore() override;

    /** This function must be implemented in each subclass to return the wavelengths and the
        corresponding relative luminosities defining the lines in the %SED. The function must
        guarantee that both arrays have the same size. The order of the pairs (other than the
        correspondance between a wavelength and its luminosity) is not important. Constant scaling
        of the luminosities is not important because the %SED will be normalized by this abstract
        class. */
    virtual void getWavelengthsAndLuminosities(Array& lambdav, Array& Lv) const = 0;

    //======================== Other Functions =======================

public:
    /** This function returns the intrinsic wavelength range of the %SED. For the current class,
        the range includes all line wavelengths. */
    Range intrinsicWavelengthRange() const override;

    /** This function returns the normalized integrated luminosity \f$L\f$ (i.e. radiative power)
        over the specified wavelength range. For the current class, this corresponds to the sum of
        the normalized luminosities for all lines with a wavelength in the specified range. */
    double integratedLuminosity(const Range& wavelengthRange) const override;

    /** This function draws a random wavelength from the normalized spectral energy distribution.
        For the current class, it chooses one of the line wavelengths. */
    double generateWavelength() const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Array _inlambdav;  // intrinsic wavelengths (i.e. as obtained from subclass)
    Array _inLv;       // intrinsic normalized specific luminosities (i.e. as obtained from subclass,
                       //                             but normalized with source range normalization)
    Array _lambdav;    // wavelengths within source range
    Array _Lv;         // normalized luminosities within source range
    Array _Pv;         // normalized cumulative distribution within source range
};

////////////////////////////////////////////////////////////////////

#endif
