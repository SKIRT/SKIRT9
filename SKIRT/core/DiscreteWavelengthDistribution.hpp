/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DISCRETEWAVELENGTHDISTRIBUTION_HPP
#define DISCRETEWAVELENGTHDISTRIBUTION_HPP

#include "WavelengthDistribution.hpp"
#include "Array.hpp"

////////////////////////////////////////////////////////////////////

/** DiscreteWavelengthDistribution is an abstract class for representing discrete wavelength
    probability distributions that are tabulated by the user in the form of a list of wavelengths.
    The wavelengths can be given in arbitrary order (they are sorted automatically) and any
    wavelengths outside of the wavelength range of the associated source (obtained through the
    WavelengthRangeInterface) are automatically removed from the list before the distribution is
    constructed. If no wavelengths are inside the source range, a fatal error is issued.

    The implemented "discrete" distribution has a constant nonzero value inside a set of narrow
    ranges placed around the tabulated wavelengths, and is zero everywhere else. The nonzero ranges
    all have the same width given by \f$w=\frac{2}{1000}\,\lambda_0\f$ where \f$\lambda_0\f$ is the
    shortest wavelength in the tabulated list. As a result, after normalization, the constant value
    of the probability distribution inside each of the ranges is given by \f$ 1 / (w\,N_\lambda)
    \f$.

    The generateWavelength() function in fact always generates exactly one of the tabulated
    wavelengths (with equal probability), rather than distributing the generated wavelengths across
    each range. The probability() function, however, returns the proper value for any given
    wavelength (i.e. nonzero within each narrow range and zero elsewehere).

    The subclass must load the tabulated wavelengths, and this abstract class handles everything
    else. */
class DiscreteWavelengthDistribution : public WavelengthDistribution
{
    ITEM_ABSTRACT(DiscreteWavelengthDistribution, WavelengthDistribution,
                  "a discrete wavelength probability distribution tabulated by the user")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function asks the subclass to load the wavelengths and precalculates some data for use
        by the other functions in this class. */
    void setupSelfBefore() override;

    /** This function must be implemented in each subclass to return the wavelengths tabulating the
        discrete distribution. The order of the wavelengths in the list does not matter. */
    virtual Array getWavelengths() const = 0;

    //======================== Other Functions =======================

public:
    /** This function returns the probability of the distribution at the given wavelength, as
        described in the class header. */
    double probability(double wavelength) const override;

    /** This function draws a random wavelength from the wavelength distribution as described in
        the class header. */
    double generateWavelength() const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Array _lambdav;             // N tabulated wavelengths
    Array _borderv;             // K=N*2 ordered range border points
    double _probability{0.};    // probability within each narrow range
};

////////////////////////////////////////////////////////////////////

#endif
