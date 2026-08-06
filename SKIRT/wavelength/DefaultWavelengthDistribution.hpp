/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DEFAULTWAVELENGTHDISTRIBUTION_HPP
#define DEFAULTWAVELENGTHDISTRIBUTION_HPP

#include "Range.hpp"
#include "WavelengthDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** A DefaultWavelengthDistribution object represents a wavelength probability distribution for
    which the logarithm of the wavelength is distributed uniformly over the full wavelength range
    of the source. This is equivalent to a LogWavelengthDistribution object configured with a
    wavelength range that fully contains the source wavelength range, but without the need for any
    user configuration.

    Specifically, with \f$a\f$ and \f$b\f$ denoting the minimum and maximum of the source
    wavelength range, the normalized probability distribution can be written as

    \f[ p(\lambda) = \begin{cases} \frac{1}{\ln b - \ln a}\,\frac{1}{\lambda} & a < \lambda < b \\
    0 & \mathrm{elsewhere} \end{cases} \f]

    Using the inversion method, this easily leads to the formula for sampling a wavelength from
    this distribution:

    \f[ \ln\lambda = \ln a + (\ln b - \ln a) \mathcal{X} \f]

    where \f$\mathcal{X} \f$ is a random deviate. */
class DefaultWavelengthDistribution : public WavelengthDistribution
{
    ITEM_CONCRETE(DefaultWavelengthDistribution, WavelengthDistribution,
                  "the default (logarithmic) wavelength probability distribution")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function precalculates some values used by the other functions in this class. */
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
    Range _range;
    double _logMin{0};
    double _logWidth{0};
};

////////////////////////////////////////////////////////////////////

#endif
