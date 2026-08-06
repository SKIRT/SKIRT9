/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef POWERLAWGRAINSIZEDISTRIBUTION_HPP
#define POWERLAWGRAINSIZEDISTRIBUTION_HPP

#include "RangeGrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** PowerLawGrainSizeDistribution is a GrainSizeDistribution subclass that represents a dust grain
    size distribution of the form \f[ \frac{\text{d}n_\text{D}}{\text{d}a} \propto a^{-\gamma}
    \qquad \text{for}\quad a_\text{min} \leq a \leq a_\text{max}, \f] where the exponent
    \f$\gamma>0\f$ can be configured as an attribute of this class, and the size range can be
    configured in the RangeGrainSizeDistribution base class.

    In fact, the function is guaranteed to exactly return the value of \f$a^{-\gamma}\f$, without
    any multiplication factor, so that the size distribution can be normalized (externally) by
    specifying a proportionality factor with the appropriate units. For example, to obtain a size
    distribution function specifying the amount of dust per hydrogen atom as a function of grain
    size, one would (externally) add a factor \f$C\f$, as in \f[
    \frac{\text{d}n_\text{D}}{\text{d}a} / n_\text{H}= C\, a^{-\gamma} \qquad \text{for}\quad
    a_\text{min} \leq a \leq a_\text{max}, \f] where \f$C\f$ has units of
    \f$\text{m}^{\gamma-1}\f$. */
class PowerLawGrainSizeDistribution : public RangeGrainSizeDistribution
{
    ITEM_CONCRETE(PowerLawGrainSizeDistribution, RangeGrainSizeDistribution, "a power-law dust grain size distribution")

        PROPERTY_DOUBLE(exponent, "the (absolute value of the) exponent in the power-law distribution function")
        ATTRIBUTE_MIN_VALUE(exponent, "]0")
        ATTRIBUTE_MAX_VALUE(exponent, "10]")
        ATTRIBUTE_DEFAULT_VALUE(exponent, "3.5")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked by classes that wish to hard-code the creation of a new
        grain size distribution object of this type (as opposed to creation through the ski file).
        Before the constructor returns, the newly created object is hooked up as a child to the
        specified parent in the simulation hierarchy (so it will automatically be deleted), its
        properties have been initialized to the specified values, and its setup() function has been
        called. */
    explicit PowerLawGrainSizeDistribution(SimulationItem* parent, double minSize, double maxSize, double exponent);

    //======================== Other Functions =======================

public:
    /** This function returns the value of \f$a^{-\gamma}\f$. */
    double dnda(double a) const override;
};

////////////////////////////////////////////////////////////////////

#endif
