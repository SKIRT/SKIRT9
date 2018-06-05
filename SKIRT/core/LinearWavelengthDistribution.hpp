/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LINEARWAVELENGTHDISTRIBUTION_HPP
#define LINEARWAVELENGTHDISTRIBUTION_HPP

#include "RangeWavelengthDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** A LinearWavelengthDistribution object represents a (uniform) linear wavelength probability
    distribution in a wavelength range configured by the user. */
class LinearWavelengthDistribution : public RangeWavelengthDistribution
{
    ITEM_CONCRETE(LinearWavelengthDistribution, RangeWavelengthDistribution,
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
