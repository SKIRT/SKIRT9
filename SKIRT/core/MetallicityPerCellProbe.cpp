/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MetallicityPerCellProbe.hpp"
#include "Configuration.hpp"
#include "MediumSystem.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"

////////////////////////////////////////////////////////////////////

void MetallicityPerCellProbe::probe()
{
    if (find<Configuration>()->hasMedium())
    {
        // locate the medium system and the grid
        auto ms = find<MediumSystem>();
        auto grid = ms->grid();
        int numMedia = ms->numMedia();
        int numCells = grid->numCells();

        // loop over the media
        for (int h = 0; h != numMedia; ++h)
        {
            // check whether this medium component has a metallicity state variable
            auto mix = ms->media()[h]->mix();
            bool hasMetallicity = false;
            for (const auto& candidate : mix->specificStateVariableInfo())
                if (candidate.identifier() == StateVariable::Identifier::Metallicity) hasMetallicity = true;

            // if there is a metallicity variable
            if (hasMetallicity)
            {
                // create a text file
                TextOutFile out(this, itemName() + "_Z_" + std::to_string(h), "metallicity");

                // write the header
                out.addColumn("spatial cell index", "", 'd');
                out.addColumn("metallicity");

                // write a line for each cell
                for (int m = 0; m != numCells; ++m)
                {
                    out.writeRow(static_cast<double>(m), ms->metallicity(m, h));
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////
