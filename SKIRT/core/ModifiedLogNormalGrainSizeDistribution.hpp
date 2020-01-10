/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MODIFIEDLOGNORMALGRAINSIZEDISTRIBUTION_HPP
#define MODIFIEDLOGNORMALGRAINSIZEDISTRIBUTION_HPP

#include "LogNormalGrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** ModifiedLogNormalGrainSizeDistribution is a GrainSizeDistribution subclass that represents a
    modified log-normal dust grain size distribution of the form \f[
    \frac{\text{d}n_\text{D}}{\text{d}a} \propto \frac{1}{a} \,\exp\left[ -
    \frac{(\ln(a/a_0))^2}{2\sigma^2} \right] \, M(a) \qquad \text{for}\quad a_\text{min} \leq a
    \leq a_\text{max} \f] with a mixing term \f[ M(a) = y_0 +
    (y_1-y_0)\frac{\ln(a/a_\text{min})}{\ln(a_\text{max}/a_\text{min})}. \f]

    The size range of the distribution can be configured in the RangeGrainSizeDistribution base
    class. The centroid \f$a_0\f$ and the width \f$\sigma\f$ can be configured in the
    LogNormalGrainSizeDistribution base class. The remaining two parameters \f$y_0\f$ and \f$y_1\f$
    can be configured as attributes in this class. The function is scaled arbitrarily.

    The functional form for the grain size distribution implemented by this class is inspired by
    the DustEM code, which is described in Compiègne et al. 2011 (AA, 525, A103) and can be
    downloaded from http://www.ias.u-psud.fr/DUSTEM/. */
class ModifiedLogNormalGrainSizeDistribution : public LogNormalGrainSizeDistribution
{
    ITEM_CONCRETE(ModifiedLogNormalGrainSizeDistribution, LogNormalGrainSizeDistribution,
                  "a modified log-normal dust grain size distribution")

        PROPERTY_DOUBLE(firstMixingParameter, "the first mixing parameter y0")
        ATTRIBUTE_MIN_VALUE(firstMixingParameter, "[0")
        ATTRIBUTE_MAX_VALUE(firstMixingParameter, "1]")
        ATTRIBUTE_DEFAULT_VALUE(firstMixingParameter, "1")

        PROPERTY_DOUBLE(secondMixingParameter, "the second mixing parameter y1")
        ATTRIBUTE_MIN_VALUE(secondMixingParameter, "[0")
        ATTRIBUTE_MAX_VALUE(secondMixingParameter, "1]")
        ATTRIBUTE_DEFAULT_VALUE(secondMixingParameter, "1")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the value of \f$\frac{\text{d}n_\text{D}}{\text{d}a}\f$ as described
        in the header for this class (with an arbitrary proportionality factor of one). */
    double dnda(double a) const override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation
    const double& _y0{_firstMixingParameter};
    const double& _y1{_secondMixingParameter};
};

////////////////////////////////////////////////////////////////////

#endif
