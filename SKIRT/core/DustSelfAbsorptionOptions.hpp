/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTSELFABSORPTIONOPTIONS_HPP
#define DUSTSELFABSORPTIONOPTIONS_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** The DustSelfAbsorptionOptions class simply offers a number of configuration options related to
    the self-consistent calculation of dust self-absorption, including the convergence criteria for
    the iteration process. */
class DustSelfAbsorptionOptions : public SimulationItem
{
    ITEM_CONCRETE(DustSelfAbsorptionOptions, SimulationItem, "a set of options related to calculating dust self-absorption")

    PROPERTY_INT(minIterations, "the minimum number of dust self-absorption iterations")
        ATTRIBUTE_MIN_VALUE(minIterations, "1")
        ATTRIBUTE_MAX_VALUE(minIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(minIterations, "1")
        ATTRIBUTE_RELEVANT_IF(minIterations, "iterateSelfAbsorption")
        ATTRIBUTE_DISPLAYED_IF(minIterations, "Level3")

    PROPERTY_INT(maxIterations, "the maximum number of dust self-absorption iterations")
        ATTRIBUTE_MIN_VALUE(maxIterations, "1")
        ATTRIBUTE_MAX_VALUE(maxIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(maxIterations, "10")
        ATTRIBUTE_RELEVANT_IF(maxIterations, "iterateSelfAbsorption")
        ATTRIBUTE_DISPLAYED_IF(maxIterations, "Level3")

    PROPERTY_DOUBLE(maxFractionOfPrimary, "convergence is reached when the total absorbed dust luminosity "
                                          "is less than this fraction of the total absorbed primary luminosity")
        ATTRIBUTE_MIN_VALUE(maxFractionOfPrimary, "]0")
        ATTRIBUTE_MAX_VALUE(maxFractionOfPrimary, "1[")
        ATTRIBUTE_DEFAULT_VALUE(maxFractionOfPrimary, "0.01")
        ATTRIBUTE_RELEVANT_IF(maxFractionOfPrimary, "iterateSelfAbsorption")
        ATTRIBUTE_DISPLAYED_IF(maxFractionOfPrimary, "Level3")

    PROPERTY_DOUBLE(maxFractionOfPrevious, "convergence is reached when the total absorbed dust luminosity "
                                           "has changed by less than this fraction compared to the previous iteration")
        ATTRIBUTE_MIN_VALUE(maxFractionOfPrevious, "]0")
        ATTRIBUTE_MAX_VALUE(maxFractionOfPrevious, "1[")
        ATTRIBUTE_DEFAULT_VALUE(maxFractionOfPrevious, "0.03")
        ATTRIBUTE_RELEVANT_IF(maxFractionOfPrevious, "iterateSelfAbsorption")
        ATTRIBUTE_DISPLAYED_IF(maxFractionOfPrevious, "Level3")

    PROPERTY_DOUBLE(iterationPacketsMultiplier,
                    "the multiplier on the number of photon packets launched for each self-absorption iteration")
        ATTRIBUTE_MIN_VALUE(iterationPacketsMultiplier, "]0")
        ATTRIBUTE_MAX_VALUE(iterationPacketsMultiplier, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(iterationPacketsMultiplier, "1")
        ATTRIBUTE_RELEVANT_IF(iterationPacketsMultiplier, "iterateSelfAbsorption")
        ATTRIBUTE_DISPLAYED_IF(iterationPacketsMultiplier, "Level3")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif
