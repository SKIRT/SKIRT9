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
    configured in the RangeGrainSizeDistribution base class. */
class PowerLawGrainSizeDistribution: public RangeGrainSizeDistribution
{
    ITEM_CONCRETE(PowerLawGrainSizeDistribution, RangeGrainSizeDistribution,
                  "a power-law dust grain size distribution")

    PROPERTY_DOUBLE(exponent, "the (absolute value of the) exponent in the power-law distribution function")
        ATTRIBUTE_MIN_VALUE(exponent, "]0")
        ATTRIBUTE_MAX_VALUE(exponent, "10]")
        ATTRIBUTE_DEFAULT_VALUE(exponent, "3.5")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the value of \f$a^{-\gamma}\f$. */
    double dnda(double a) const override;
};

////////////////////////////////////////////////////////////////////

#endif
