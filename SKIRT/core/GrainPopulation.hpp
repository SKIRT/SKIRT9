/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GRAINPOPULATION_HPP
#define GRAINPOPULATION_HPP

#include "GrainComposition.hpp"
#include "GrainSizeDistribution.hpp"
#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** GrainPopulation is simple class that represents a particular dust grain population. A grain
    population is essentially defined by the combination of a grain composition (an instance of a
    GrainComposition subclass) providing the optical and calorimetric properties of the grain
    material, and a grain size distribution (an instance of a GrainSizeDistribution subclass). In
    addition, the amount of dust contained in the population is specified, in one of three ways: as
    an absolute dust mass per hydrogen atom, as a ratio of dust mass per hydrogen mass, or by using
    a given proportionality factor on the size distribution (which in that case should have some
    known normalization).

    Finally, a grain population specifies the number of grain size bins (on a logarithmic scale)
    that should be used when calculating emission spectra for the grain population. Because the
    emissivity is nonlinear as a function of grain size, the calculation is performed for a
    representative grain in each bin. A larger number of bins results in higher accuracy. On the
    other hand, the required computation time and memory consumption increase roughly linearly with
    the number of bins. */
class GrainPopulation : public SimulationItem
{
    /** The enumeration type indicating the mechanism for specifying the amount of dust in the
        population. */
    ENUM_DEF(NormalizationType, DustMassPerHydrogenAtom, DustMassPerHydrogenMass, FactorOnSizeDistribution)
        ENUM_VAL(NormalizationType, DustMassPerHydrogenAtom, "an absolute dust mass per hydrogen atom")
        ENUM_VAL(NormalizationType, DustMassPerHydrogenMass, "a ratio of dust mass per hydrogen mass")
        ENUM_VAL(NormalizationType, FactorOnSizeDistribution, "a proportionality factor on the size distribution")
    ENUM_END()

    ITEM_CONCRETE(GrainPopulation, SimulationItem, "a dust grain population")

        PROPERTY_ITEM(composition, GrainComposition, "the dust grain composition")
        ATTRIBUTE_DEFAULT_VALUE(composition, "DraineGraphiteGrainComposition")

        PROPERTY_ITEM(sizeDistribution, GrainSizeDistribution, "the dust grain size distribution")
        ATTRIBUTE_DEFAULT_VALUE(sizeDistribution, "PowerLawGrainSizeDistribution")

        PROPERTY_INT(numSizes, "the number of grain size bins")
        ATTRIBUTE_MIN_VALUE(numSizes, "1")
        ATTRIBUTE_MAX_VALUE(numSizes, "1000")
        ATTRIBUTE_DEFAULT_VALUE(numSizes, "8")

        PROPERTY_ENUM(normalizationType, NormalizationType,
                      "the mechanism for specifying the amount of dust in the population")
        ATTRIBUTE_DEFAULT_VALUE(normalizationType, "DustMassPerHydrogenMass")

        PROPERTY_DOUBLE(dustMassPerHydrogenAtom, "the dust mass per hydrogen atom")
        ATTRIBUTE_QUANTITY(dustMassPerHydrogenAtom, "mass")
        ATTRIBUTE_MIN_VALUE(dustMassPerHydrogenAtom, "[0")
        ATTRIBUTE_RELEVANT_IF(dustMassPerHydrogenAtom, "normalizationTypeDustMassPerHydrogenAtom")

        PROPERTY_DOUBLE(dustMassPerHydrogenMass, "the dust mass per hydrogen mass")
        ATTRIBUTE_MIN_VALUE(dustMassPerHydrogenMass, "[0")
        ATTRIBUTE_MAX_VALUE(dustMassPerHydrogenMass, "1[")
        ATTRIBUTE_RELEVANT_IF(dustMassPerHydrogenMass, "normalizationTypeDustMassPerHydrogenMass")

        PROPERTY_DOUBLE(factorOnSizeDistribution, "the proportionality factor on the size distribution")
        ATTRIBUTE_MIN_VALUE(factorOnSizeDistribution, "[0")
        ATTRIBUTE_DEFAULT_VALUE(factorOnSizeDistribution, "1")
        ATTRIBUTE_RELEVANT_IF(factorOnSizeDistribution, "normalizationTypeFactorOnSizeDistribution")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked by classes that wish to hard-code the creation of a new
        grain population object (as opposed to creation through the ski file). Before the
        constructor returns, the newly created object is hooked up as a child to the specified
        parent in the simulation hierarchy (so it will automatically be deleted), its properties
        have been set to the values specified in the constructor, and its setup() function has been
        called. The caller must guarantee that the lifetime of the specified composition and size
        distribution objects is as least as long as the lifetime of the newly created grain
        population object.

        The \em normValue argument specifies the normalization value corresponding to the specified
        normalization type, i.e. dustMassPerHydrogenAtom (in kg), dustMassPerHydrogenMass
        (dimensionless ratio), or factorOnSizeDistribution (dimensionless factor). */
    explicit GrainPopulation(SimulationItem* parent, GrainComposition* composition,
                             GrainSizeDistribution* sizeDistribution, int numSizes, NormalizationType normType,
                             double normValue);
};

////////////////////////////////////////////////////////////////////

#endif
