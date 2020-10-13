/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustDestructionRecipe.hpp"
#include "Configuration.hpp"
#include "FragmentDustMixDecorator.hpp"
#include "Log.hpp"
#include "MaterialState.hpp"
#include "MediumSystem.hpp"

////////////////////////////////////////////////////////////////////

void DustDestructionRecipe::beginUpdate(int /*numCells*/)
{
    // store radiation field wavelength grid
    _dlambdav = find<Configuration>()->radiationFieldWLG()->dlambdav();

    // store a pointer for each medium component: either the corresponding fragmented dust mix or a null pointer
    for (auto medium : find<MediumSystem>()->media())
    {
        auto fd = dynamic_cast<const FragmentDustMixDecorator*>(medium->mix());
        if (fd && !fd->hasDynamicDensities())
        {
            find<Log>()->warning("Cannot destroy dust for fragmented dust mix without dynamic densities");
            fd = nullptr;
        }
        _fdv.push_back(fd);
    }
}

////////////////////////////////////////////////////////////////////

bool DustDestructionRecipe::update(MaterialState* state, const Array& Jv)
{
    // flag becomes true whenever the state is updated
    bool updated = false;

    // act only if this component has a fragmented dust mix and if this cell contains dust
    auto fd = _fdv[state->mediumIndex()];
    if (fd && state->numberDensity() > 0.)
    {
        // loop over all fragments
        int numFrags = fd->numPopulations();
        for (int f = 0; f != numFrags; ++f)
        {
            // act only if this fragment contains dust
            if (state->custom(f) > 0.)
            {
                // get the grain type and representative grain mass for this fragment
                bool graphite = fd->populationIsGraphite(f);
                double m = fd->populationGrainMass(f);

                // get the equilibrium temperature for this fragment in the cell's radiation field
                double T = fd->populationTemperature(f, Jv);

                // call the destruction recipe to determine the density fraction
                double fraction = densityFraction(graphite, m, T);

                // if the new fraction differs substantially from the previous one, update the medium state
                if (abs(fraction - state->custom(numFrags + f)) > 0.01)  // TO DO: allow configuring convergence
                {
                    state->setCustom(numFrags + f, fraction);
                    updated = true;
                }
            }
        }
    }
    return updated;
}

////////////////////////////////////////////////////////////////////

double DustDestructionRecipe::densityFraction(bool graphite, double mass, double temperature)
{
    return 0.5;
}

////////////////////////////////////////////////////////////////////
