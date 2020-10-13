/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTDESTRUCTIONRECIPE_HPP
#define DUSTDESTRUCTIONRECIPE_HPP

#include "DynamicStateRecipe.hpp"
class FragmentDustMixDecorator;

////////////////////////////////////////////////////////////////////

/** DustDestructionRecipe is a dynamic medium state recipe TO DO ... */
class DustDestructionRecipe : public DynamicStateRecipe
{
    ITEM_CONCRETE(DustDestructionRecipe, DynamicStateRecipe, "a dust destruction recipe")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function is called before an update cycle begins. It caches the radiation field
        wavelength grid and information on the dust population fragments in the medium system
        offered by one or more FragmentDustMixDecorator instances. */
    void beginUpdate(int numCells) override;

    /** This function is called repeatedly as part of the update cycle. If the density of the
        medium component in the cell under study is nonzero, it calculates and if necessary updates
        the dynamic density fractions in the corresponding medium state. TO DO ...
        The function returns true if the density has been updated and false otherwise. */
    bool update(MaterialState* state, const Array& Jv) override;

    /** This function returns the density fraction for a grain population with the specified type
        (graphite or silicate), representative grain mass and equilibrium temperature. TO DO ... */
    virtual double densityFraction(bool graphite, double mass, double temperature);

    //======================== Data Members =======================

private:
    Array _dlambdav;                               // wavelength bin widths for the radiation field wavelength grid
    vector<const FragmentDustMixDecorator*> _fdv;  // fragmented mixes in the simulation (or nullptr for other mix)
};

////////////////////////////////////////////////////////////////////

#endif
