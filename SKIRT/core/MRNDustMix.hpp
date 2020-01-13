/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MRNDUSTMIX_HPP
#define MRNDUSTMIX_HPP

#include "MultiGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The MRNDustMix class represents a dust mixture consisting of silicate graphite and dust grains
    with optical properties taken from Bruce Draine's web site and with grain size distributions
    according to the well-known MRN power-law prescription (Mathis, Rumpl & Nordsieck 1977, ApJ,
    217, 425). The actual coefficients for the size distributions are taken from Weingartner &
    Draine (2001, ApJ, 548, 296) on page 296, i.e. \f[ \frac{\text{d}n_\text{D}}{\text{d}a} /
    n_\text{H}= C\, a^{-3.5} \qquad \text{for}\quad a_\text{min} \leq a \leq a_\text{max}, \f] with
    \f$C=10^{-25.13}\,\text{cm}^{2.5}\f$ for graphite and \f$C=10^{-25.11}\,\text{cm}^{2.5}\f$ for
    silicate, and with \f$a_\text{min}=50\,\text{\AA}\f$ and \f$a_\text{max}=0.25\,\mu\text{m}\f$
    for both grain populations.

    For each of the populations, the user can configure the number of grain size bins used to
    discretize the thermal emission calculations for the dust mix.

    The class uses the grain compositions defined in the DraineSilicateGrainComposition and
    DraineGraphiteGrainComposition classes, and the power-law size distribution defined in the
    PowerLawGrainSizeDistribution class. */
class MRNDustMix : public MultiGrainDustMix
{
    ITEM_CONCRETE(MRNDustMix, MultiGrainDustMix, "an MRN (1997) dust mix")

        PROPERTY_INT(numSilicateSizes, "the number of silicate grain size bins")
        ATTRIBUTE_MIN_VALUE(numSilicateSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numSilicateSizes, "5")

        PROPERTY_INT(numGraphiteSizes, "the number of graphite grain size bins")
        ATTRIBUTE_MIN_VALUE(numGraphiteSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numGraphiteSizes, "5")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function adds the relevant grain populations to the dust mix */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
