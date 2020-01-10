/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LOGNORMALGRAINSIZEDISTRIBUTION_HPP
#define LOGNORMALGRAINSIZEDISTRIBUTION_HPP

#include "RangeGrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** LogNormalGrainSizeDistribution is a GrainSizeDistribution subclass that represents a log-normal
    dust grain size distribution of the form \f[ \frac{\text{d}n_\text{D}}{\text{d}a} \propto
    \frac{1}{a} \,\exp\left[ - \frac{(\ln(a/a_0))^2}{2\sigma^2} \right] \qquad \text{for}\quad
    a_\text{min} \leq a \leq a_\text{max}. \f]

    The size range of the distribution can be configured in the RangeGrainSizeDistribution base
    class. The remaining two parameters, the centroid \f$a_0\f$ and the width \f$\sigma\f$, can be
    configured as attributes in this class. The function is scaled arbitrarily.

    The functional form for the grain size distribution implemented by this class is inspired by
    the DustEM code, which is described in Compiègne et al. 2011 (AA, 525, A103) and can be
    downloaded from http://www.ias.u-psud.fr/DUSTEM/. */
class LogNormalGrainSizeDistribution : public RangeGrainSizeDistribution
{
    ITEM_CONCRETE(LogNormalGrainSizeDistribution, RangeGrainSizeDistribution,
                  "a log-normal dust grain size distribution")

        PROPERTY_DOUBLE(centroid, "the centroid a0 of the log-normal law")
        ATTRIBUTE_QUANTITY(centroid, "grainsize")
        ATTRIBUTE_MIN_VALUE(centroid, "]0")
        ATTRIBUTE_MAX_VALUE(centroid, "1 mm]")
        ATTRIBUTE_DEFAULT_VALUE(centroid, "1 nm")

        PROPERTY_DOUBLE(width, "the width σ of the log-normal law")
        ATTRIBUTE_MIN_VALUE(width, "]0")
        ATTRIBUTE_MAX_VALUE(width, "5]")
        ATTRIBUTE_DEFAULT_VALUE(width, "0.4")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked by classes that wish to hard-code the creation of a new
        grain size distribution object of this type (as opposed to creation through the ski file).
        Before the constructor returns, the newly created object is hooked up as a child to the
        specified parent in the simulation hierarchy (so it will automatically be deleted), its
        properties have been initialized to the specified values, and its setup() function has been
        called. */
    explicit LogNormalGrainSizeDistribution(SimulationItem* parent, double minSize, double maxSize, double centroid,
                                            double width);

    //======================== Other Functions =======================

public:
    /** This function returns the value of \f$\frac{\text{d}n_\text{D}}{\text{d}a}\f$ as described
        in the header for this class (with an arbitrary proportionality factor of one). */
    double dnda(double a) const override;
};

////////////////////////////////////////////////////////////////////

#endif
