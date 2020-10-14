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
    //============= Identity =============

protected:
    /** This function returns the name of the item's type and is usually automatically provided
        through the ITEM macro's that define the metadata for a simulation item. Because
        DustMixFragment cannot be configured in a ski file, that metadata is not needed and we
        explicitly define the type() function instead. It returns the name of the
        FragmentDustMixDecorator class for which DustMixFragment serves under the hood. */
    string type() const override;

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
        the specified single grain population, limited to the specified grain size range. The last
        argument specifies the grain size distribution normalization factor for the population. The
        GrainPopulation object must be kept alive by the caller as long as the DustMixFragment
        instance lives.

        Before the constructor returns, the newly created object is hooked up as a child to the
        specified parent in the simulation hierarchy (so it will automatically be deleted), the
        size-limited grain population has been passed to the MultiGrainDustMix base class, and the
        setup() function has been called. */
    explicit DustMixFragment(SimulationItem* parent, ScatteringMode scatteringMode, const GrainPopulation* population,
                             double amin, double amax, double normalization);

protected:
    /** This function adds the grain population specified in the constructor to the multi-grain
        dust mix through the appropriate base class function. */
    void setupSelfBefore() override;

    /** This function calculates and stores some derived properties of the dust population
        represented by this fragment, intended for use by dust destruction/allocation recipes. */
    void setupSelfAfter() override;

    //======== Scattering implementation for dust mixes =======

public:
    /** This function returns the scattering mode supported by this dust mix fragment as passed to
        the constructor. */
    ScatteringMode scatteringMode() const override;

    //============= Extra queries for cached information =============

public:
    /** This function returns true if the human-readable identifier for the type of grain material
        represented by this fragment contains "Gra", "PAH" or "CM20", and false otherwise. */
    bool isGraphite() const;

    /** This function returns the average radius of a dust grain in the population represented by
        this fragment. */
    double grainRadius() const;

    //======================== Data Members ========================

private:
    // initialized by constructor
    ScatteringMode _scatteringMode{ScatteringMode::HenyeyGreenstein};
    const GrainPopulation* _population{nullptr};

    // initialized by setupSelfAfter
    bool _isGraphite{false};
    double _grainRadius{0.};
};

////////////////////////////////////////////////////////////////////

#endif
