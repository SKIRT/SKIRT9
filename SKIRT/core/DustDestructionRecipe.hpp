/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTDESTRUCTIONRECIPE_HPP
#define DUSTDESTRUCTIONRECIPE_HPP

#include "DynamicStateRecipe.hpp"
class FragmentDustMixDecorator;

////////////////////////////////////////////////////////////////////

/** DustDestructionRecipe is an abstract class representing recipes for dynamically calculating
    radiative destruction of dust grains by iterating over the radiation field. The class derives
    from DynamicStateRecipe so that it fits in the overall dynamic medium state framework. It
    cooperates with the FragmentDustMixDecorator class to handle each dust mix fragment separately
    (i.e. different grain materials and/or grain size bins). A subclass can specify a concrete
    destruction recipe by implementing a single function that takes properties of the radiation
    field and grain population/size bin under study and returns the corresponding non-destroyed
    density fraction.

    When a DustDestructionRecipe instance is included in the dynamic state configuration of a
    simulation, it requires there to be one or more medium components with a
    FragmentDustMixDecorator material mix that has the \em hasDynamicDensities flag enabled. These
    media components will be automatically detected and handled. */
class DustDestructionRecipe : public DynamicStateRecipe
{
    ITEM_ABSTRACT(DustDestructionRecipe, DynamicStateRecipe, "a dust destruction recipe")

        PROPERTY_DOUBLE(densityFractionTolerance, "the convergence tolerance on the dynamic density fraction")
        ATTRIBUTE_MIN_VALUE(densityFractionTolerance, "[1e-6")
        ATTRIBUTE_MAX_VALUE(densityFractionTolerance, "0.5]")
        ATTRIBUTE_DEFAULT_VALUE(densityFractionTolerance, "0.05")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function is called before an update cycle begins. It caches the radiation field
        wavelength grid and information on the dust population fragments in the medium system
        offered by one or more FragmentDustMixDecorator instances. */
    void beginUpdate(int numCells) override;

    /** This function is called repeatedly as part of the update cycle. If the density of the
        medium component in the cell under study is nonzero, it calculates the non-destroyed
        density fraction for each of the dust mix fragments in the medium component under study by
        calling the densityFraction() function defined in a subclass. If a newly calculated
        fraction differs by more than some small tolerance from the previously stored value, the
        corresponding medium state is updated. The function returns \em NotUpdated if none of the
        fractions have been updated, \em UpdatedNotConverged if one or more fractions were updated
        by more than the configured convergence tolerance, \em UpdatedConverged if one or more
        fractions were updated but none of them differed by more than the configured convergence
        tolerance. */
    UpdateStatus update(MaterialState* state, const Array& Jv) override;

    /** This function returns the non-destroyed density fraction for a grain population with the
        specified type (graphite or silicate), average grain radius, radiation field, and
        equilibrium temperature. It must be implemented by a subclass. */
    virtual double densityFraction(bool graphite, double a, const Array& Jv, double T) const = 0;

    //======================== Data Members =======================

private:
    Array _dlambdav;                               // wavelength bin widths for the radiation field wavelength grid
    vector<const FragmentDustMixDecorator*> _fdv;  // fragmented mixes in the simulation (or nullptr for other mix)
};

////////////////////////////////////////////////////////////////////

#endif
