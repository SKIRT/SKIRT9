/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTMIXFRAGMENT_HPP
#define DUSTMIXFRAGMENT_HPP

#include "MultiGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** DustMixFragment is a helper class that represents a single grain population created by another
    dust mix. It is used by the FragmentDustMixDecorator class to represent a \em fragment of a
    MultiGrainDustMix instance.

    DustMixFragment itself inherits from MultiGrainDustMix as well, so that it can use the full
    power of this base class (albeit for just a single dust population). Although DustMixFragment
    indirectly inherits from SimulationItem, and thus can be part of the simulation hierarchy, its
    instances can only be constructed programmatically as opposed to being configured as part of a
    ski file. */
class DustMixFragment : public MultiGrainDustMix
{
    //============= Construction - Setup - Destruction =============

public:
    /** This constructor programmatically constructs a DustMixFragment instance from an existing
        grain population (there are no provisions for creation through the ski file). The resulting
        object represents a fully functional "multi-grain" dust mix defined by the specified single
        grain population. The GrainPopulation object must be kept alive by the caller as long as
        the DustMixFragment instance lives.

        Before the constructor returns, the newly created object is hooked up as a child to the
        specified parent in the simulation hierarchy (so it will automatically be deleted), the
        specified grain population has been passed to the MultiGrainDustMix base class, and the
        setup() function has been called. */
    explicit DustMixFragment(SimulationItem* parent, ScatteringMode scatteringMode, const GrainPopulation* population);

    /** This constructor programmatically constructs a DustMixFragment instance from a given size
        bin of an existing grain population (there are no provisions for creation through the ski
        file). The resulting object represents a fully functional "multi-grain" dust mix defined by
        the specified single grain population, limited to the specified grain size range. The
        GrainPopulation object must be kept alive by the caller as long as the DustMixFragment
        instance lives.

        Before the constructor returns, the newly created object is hooked up as a child to the
        specified parent in the simulation hierarchy (so it will automatically be deleted), the
        size-limited grain population has been passed to the MultiGrainDustMix base class, and the
        setup() function has been called. */
    explicit DustMixFragment(SimulationItem* parent, ScatteringMode scatteringMode, const GrainPopulation* population,
                             double amin, double amax);

protected:
    /** This function adds the grain population specified in the constructor to the multi-grain
        dust mix through the appropriate base class function. */
    void setupSelfBefore() override;

    //======== Scattering implementation for dust mixes =======

public:
    /** This function returns the scattering mode supported by this dust mix fragment as passed to
        the constructor. */
    ScatteringMode scatteringMode() const override;

    //======================== Data Members ========================

private:
    ScatteringMode _scatteringMode{ScatteringMode::HenyeyGreenstein};
    const GrainPopulation* _population{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
