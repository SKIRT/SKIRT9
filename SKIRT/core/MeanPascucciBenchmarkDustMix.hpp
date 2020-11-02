/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEANPASCUCCIBENCHMARKDUSTMIX_HPP
#define MEANPASCUCCIBENCHMARKDUSTMIX_HPP

#include "SingleGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The MeanPascucciBenchmarkDustMix class represents a population of identical dust grains used
    for the 2D radiative transfer benchmark calculations of Pascucci et al. (2004, A&A, 417, 793).
    It consists of spherical astronomical silicate grains with a grain size of 0.12 micron.
    Scattering is assumed to be isotropic. The extinction coefficients in the benchmark data are
    scale-free; we arbitrarily scale them to a reasonable order of magnitude.

    The optSi.dat data file has been downloaded from
    http://www.mpia.de/PSF/PSFpages/RT/benchmark.html at MPIA, where also additional information on
    the 2D benchmark models can be found. */
class MeanPascucciBenchmarkDustMix : public SingleGrainDustMix
{
    ITEM_CONCRETE(MeanPascucciBenchmarkDustMix, SingleGrainDustMix,
                  "a Pascucci 2D benchmark dust mix (mean properties)")
        ATTRIBUTE_TYPE_DISPLAYED_IF(MeanPascucciBenchmarkDustMix, "Level2")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function returns the name of the stored table resource tabulating the basic optical
        properties for this dust mix. */
    string resourceNameForOpticalProps() const override;
};

////////////////////////////////////////////////////////////////////

#endif
