/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SINGLEGRAINSIZEDISTRIBUTION_HPP
#define SINGLEGRAINSIZEDISTRIBUTION_HPP

#include "GrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** SingleGrainSizeDistribution represents a grain size distribution approximating a delta function
    at some specific grain size. The single grain size is configured through an attribute managed
    by this class. The amin() and amax() functions in this class return a very narrow range of
    width \f$\Delta a=a_\text{s}/1000\f$ centered on the specified size \f$a_\text{s}\f$, and the
    function dnda() returns a constant distribution function value (arbitrarily scaled). */
class SingleGrainSizeDistribution : public GrainSizeDistribution
{
    ITEM_CONCRETE(SingleGrainSizeDistribution, GrainSizeDistribution, "a single-size dust grain size distribution")

        PROPERTY_DOUBLE(size, "the single grain size for this distribution")
        ATTRIBUTE_QUANTITY(size, "grainsize")
        ATTRIBUTE_MIN_VALUE(size, "[1 Angstrom")
        ATTRIBUTE_MAX_VALUE(size, "1 mm]")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the minimum grain size \f$a_\text{min} = a_\text{s} -
        \frac{1}{2}\Delta a\f$ with \f$\Delta a=a_\text{s}/1000\f$. */
    double amin() const override;

    /** This function returns the maximum grain size \f$a_\text{max} = a_\text{s} +
        \frac{1}{2}\Delta a\f$ with \f$\Delta a=a_\text{s}/1000\f$. */
    double amax() const override;

    /** This function returns the constant value distribution function value, arbitrarily set to
        one. */
    double dnda(double a) const override;
};

////////////////////////////////////////////////////////////////////

#endif
