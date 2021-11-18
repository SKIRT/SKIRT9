/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TABULATEDWAVELENGTHDISTRIBUTION_HPP
#define TABULATEDWAVELENGTHDISTRIBUTION_HPP

#include "Array.hpp"
#include "WavelengthDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** TabulatedWavelengthDistribution is an abstract class for representing wavelength probability
    distributions that are tabulated by the user in the form of wavelength/probability pairs. The
    probability distribution function is defined segment-wise by the tabulated values, using
    logarithmic interpolation. The wavelengths must be listed in increasing or decreasing order.
    Probability values outside the range indicated by the first and the last tabulated wavelength
    are considered to be zero. In addition, this range is intersected with the wavelength range of
    the associated source (obtained through the SourceWavelengthRangeInterface) before the
    distribution is normalized.

    The subclass must load the tabulated data, and this abstract class handles everything else. */
class TabulatedWavelengthDistribution : public WavelengthDistribution
{
    ITEM_ABSTRACT(TabulatedWavelengthDistribution, WavelengthDistribution,
                  "a wavelength probability distribution tabulated by the user")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function asks the subclass to load the wavelength/probability pairs and precalculates
        the cumulative distribution for use by the other functions in this class. */
    void setupSelfBefore() override;

    /** This function must be implemented in each subclass to return the wavelengths and the
        corresponding probabilities tabulating the distribution. The function must guarantee that
        both arrays have the same size. The wavelengths must be listed in increasing or decreasing
        order. Constant scaling of the probabilities is not important because the distribution will
        be normalized by this abstract class. */
    virtual void getWavelengthsAndProbabilities(Array& lambdav, Array& pv) const = 0;

    //======================== Other Functions =======================

public:
    /** This function returns the probability of the distribution at the given wavelength. */
    double probability(double wavelength) const override;

    /** This function draws a random wavelength from the wavelength distribution. */
    double generateWavelength() const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Array _lambdav;  // wavelengths
    Array _pv;       // probability distribution, normalized to unity
    Array _Pv;       // cumulative probability distribution
};

////////////////////////////////////////////////////////////////////

#endif
