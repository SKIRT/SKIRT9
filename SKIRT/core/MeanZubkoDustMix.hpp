/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEANZUBKODUSTMIX_HPP
#define MEANZUBKODUSTMIX_HPP

#include "SingleGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The MeanZubkoDustMix class represents a population of identical dust grains with properties
    approximating those of a mixture of bare (i.e. non-composite) graphite, silicate and PAH dust
    grains. The size distribution of each of these dust grain populations is fine-tuned in such a
    way that the global dust properties accurately reproduce the extinction, emission and abundance
    constraints on the Milky Way. The size distributions are taken from Zubko, Dwek & Arendt (2004,
    ApJS, 152, 211) and correspond to model BARE_GR_S.

    The dust mass per hydrogen nucleon is obtained from table 6 of the Zubko et al. 2004 paper. The
    source of the data in electronic form is unkown. */
class MeanZubkoDustMix : public SingleGrainDustMix
{
    ITEM_CONCRETE(MeanZubkoDustMix, SingleGrainDustMix, "a Zubko et al. 2004 dust mix (mean properties)")
    ITEM_END()

    //======================== Other Functions =======================

    /** This function returns the scattering mode supported by this material mix. For this dust
        mix, it always returns the basic HenyeyGreenstein scattering mode. */
    ScatteringMode scatteringMode() const override;

    /** This function returns the name of the stored table resource tabulating the basic optical
        properties for this dust mix. */
    string resourceNameForOpticalProps() const override;
};

////////////////////////////////////////////////////////////////////

#endif
