/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTGRAINSIZEDISTRIBUTIONPROBE_HPP
#define DUSTGRAINSIZEDISTRIBUTIONPROBE_HPP

#include "SpecialtyProbe.hpp"

////////////////////////////////////////////////////////////////////

/** DustGrainSizeDistributionProbe outputs column text files tabulating the dust grain size
    distributions of the applicable dust grain populations configured in the simulation. For each
    medium component, the probe retrieves a representative material mix (the mix at the origin of
    the model coordinate system). If the material mix offers MultiGrainDustMix capabilities, i.e.
    if it is described by one or more individual dust grain populations, the probe creates a size
    distribution file for each of those populations. The files are named
    <tt>prefix_probe_grainsizes_N_M.fits</tt> where N is replaced with the zero-based index of the
    medium in the configuration and M is replaced with the zero-based index of the grain population
    in the corresponding material mix.

    Each text file tabulates \f$\text{dnda}(a)\f$ (second column) as a function of \f$a\f$ (first
    column). In other words, it tabulates the relative number of dust grains with size \f$a\f$ in
    the population, \f[ \text{dnda}(a) \propto \frac{\text{d}n_\text{D}}{\text{d}a} \qquad
    \text{for}\quad a_\text{min} \leq a \leq a_\text{max}. \f] The probe reports the
    \f$\text{dnda}(a)\f$ values as returned by the grain size distribution, without extra scaling
    or normalization. The grain size values in the table are distributed logarithmically over the
    range of the size distribution. The number of values in the table can be configured by the
    user. */
class DustGrainSizeDistributionProbe : public SpecialtyProbe
{
    ITEM_CONCRETE(DustGrainSizeDistributionProbe, SpecialtyProbe, "properties: dust grain size distribution")
        ATTRIBUTE_TYPE_DISPLAYED_IF(DustGrainSizeDistributionProbe, "Level2&Medium&MultiGrainDustMix")

        PROPERTY_INT(numSamples, "the number of samples in the size distribution table")
        ATTRIBUTE_MIN_VALUE(numSamples, "3")
        ATTRIBUTE_MAX_VALUE(numSamples, "100000")
        ATTRIBUTE_DEFAULT_VALUE(numSamples, "250")

    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif
