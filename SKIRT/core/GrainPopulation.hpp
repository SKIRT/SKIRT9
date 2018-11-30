/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GRAINPOPULATION_HPP
#define GRAINPOPULATION_HPP

#include "SimulationItem.hpp"
#include "GrainComposition.hpp"
#include "GrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** GrainPopulation is simple class that represents a particular dust grain population. A grain
    population is essentially defined by the combination of a grain composition (an instance of a
    GrainComposition subclass) providing the optical and calorimetric properties of the grain
    material, and a grain size distribution (an instance of a GrainSizeDistribution subclass). In
    addition, a grain population instance provides a relative mass weight factor. This factor is
    used when multiple grain populations are combined in a dust mix to determine the fraction of
    the total dust mass per hydrogen atom represented by each grain population.

    Finally, a grain population specifies the number of grain size bins (on a logarithmic scale)
    that should be used when calculating emission spectra for the grain population. Because the
    emissivity is nonlinear as a function of grain size, the calculation is performed for a
    representative grain in each bin. A larger number of bins results in higher accuracy. On the
    other hand, the required computation time and memory consumption increase roughly linearly with
    the number of bins. */
class GrainPopulation : public SimulationItem
{
    ITEM_CONCRETE(GrainPopulation, SimulationItem, "a dust grain population")

    PROPERTY_ITEM(composition, GrainComposition, "the dust grain composition")
        ATTRIBUTE_DEFAULT_VALUE(composition, "DraineGraphiteGrainComposition")

    PROPERTY_ITEM(sizeDistribution, GrainSizeDistribution, "the dust grain size distribution")
        ATTRIBUTE_DEFAULT_VALUE(sizeDistribution, "PowerLawGrainSizeDistribution")

    PROPERTY_DOUBLE(massWeight, "the mass weight of this grain population in the dust mix")
        ATTRIBUTE_MIN_VALUE(massWeight, "]0")
        ATTRIBUTE_MAX_VALUE(massWeight, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(massWeight, "1")

    PROPERTY_INT(numSizebins, "the number of grain size bins")
        ATTRIBUTE_MIN_VALUE(numSizebins, "1")
        ATTRIBUTE_MAX_VALUE(numSizebins, "1000")
        ATTRIBUTE_DEFAULT_VALUE(numSizebins, "8")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif
