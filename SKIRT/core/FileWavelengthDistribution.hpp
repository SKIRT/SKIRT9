/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILEWAVELENGTHDISTRIBUTION_HPP
#define FILEWAVELENGTHDISTRIBUTION_HPP

#include "TabulatedWavelengthDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** A FileWavelengthDistribution object represents a wavelength probability distribution that is
    loaded from an input file. The floating point numbers in the first two columns of the text file
    specify respectively the wavelength and the corresponding probability value. Any additional
    columns in the file are ignored. The probability distribution function is defined segment-wise
    by the tabulated values, using logarithmic interpolation. The wavelengths must be listed in
    increasing or decreasing order. Probability values outside the range indicated by the first and
    the last wavelength in the file are considered to be zero.

    The wavelengths are by default given in micron (the units can be overridden by column header
    info in the file). The probability values are in fact given in luminosity units. The default is
    to use per-wavelength units. Again, this can be overridden by column header info in the file to
    use per-frequency or neutral units. Other than this, the scaling of the values is arbitrary
    because the distribution will be normalized after being loaded. However, the input procedure
    still insists on knowing the precise units. */
class FileWavelengthDistribution : public TabulatedWavelengthDistribution
{
    ITEM_CONCRETE(FileWavelengthDistribution, TabulatedWavelengthDistribution,
                  "a wavelength probability distribution loaded from a text file")

        PROPERTY_STRING(filename, "the name of the file with the wavelength probability distribution")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function loads the input file into the specified arrays. */
    void getWavelengthsAndProbabilities(Array& lambdav, Array& pv) const override;
};

////////////////////////////////////////////////////////////////////

#endif
