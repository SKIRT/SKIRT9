/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LINWAVELENGTHDISTRIBUTION_HPP
#define LINWAVELENGTHDISTRIBUTION_HPP

#include "RangeWavelengthDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** A LinWavelengthDistribution object represents a (uniform) linear wavelength probability
    distribution in a wavelength range configured by the user. */
class LinWavelengthDistribution : public RangeWavelengthDistribution
{
    ITEM_CONCRETE(LinWavelengthDistribution, RangeWavelengthDistribution,
                  "a linear wavelength probability distribution")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the probability of the distribution at the given wavelength. */
    double probability(double wavelength) const override;

    /** This function draws a random wavelength from the wavelength distribution. */
    double generateWavelength() const override;
};

////////////////////////////////////////////////////////////////////

#endif
