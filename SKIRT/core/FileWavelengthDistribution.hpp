/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILEWAVELENGTHDISTRIBUTION_HPP
#define FILEWAVELENGTHDISTRIBUTION_HPP

#include "WavelengthDistribution.hpp"
#include "Array.hpp"

////////////////////////////////////////////////////////////////////

/** A FileWavelengthDistribution object represents a wavelength probability distribution that is
    loaded from an input file. The floating point numbers in the first two columns of the text file
    specify respectively the wavelength and the corresponding probability value. Any additional
    columns in the file are ignored. The wavelengths must be given in micron and must listed be in
    increasing order. The probability values can be given in arbitrary units because the
    distrubition will be normalized after being loaded. Probability values outside the
    range indicated by the first and the last wavelength in the file are considered to be zero. */
class FileWavelengthDistribution : public WavelengthDistribution
{
    ITEM_CONCRETE(FileWavelengthDistribution, WavelengthDistribution,
                  "a wavelength probability distribution loaded from a text file")

    PROPERTY_STRING(filename, "the name of the file with the wavelength probability distribution")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function loads the input file and precalculates the cumulative distribition for use by
        the other functions in this class. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the probability of the distribution at the given wavelength. */
    double probability(double wavelength) const override;

    /** This function draws a random wavelength from the wavelength distribution. */
    double generateWavelength() const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Array _lambdav;     // wavelengths
    Array _pv;          // probability distribution, normalized to unity
    Array _Pv;          // cumulative probability distribution
};

////////////////////////////////////////////////////////////////////

#endif
