/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CustomStatePerCellProbe.hpp"
#include "Configuration.hpp"
#include "MediumSystem.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void CustomStatePerCellProbe::probeSetup()
{
    if (probeAfter() == ProbeAfter::Setup) probe();
}

////////////////////////////////////////////////////////////////////

void CustomStatePerCellProbe::probeRun()
{
    if (probeAfter() == ProbeAfter::Run) probe();
}

////////////////////////////////////////////////////////////////////

void CustomStatePerCellProbe::probe()
{
    if (find<Configuration>()->hasMedium())
    {
        // locate the medium system and the grid
        auto units = find<Units>();
        auto ms = find<MediumSystem>();
        auto grid = ms->grid();
        int numMedia = ms->numMedia();
        int numCells = grid->numCells();

        // loop over the media
        for (int h = 0; h != numMedia; ++h)
        {
            // get the custom state variable descriptors for this medium
            auto mix = ms->media()[h]->mix();
            vector<StateVariable> descriptors;
            for (auto candidate : mix->specificStateVariableInfo())
                if (candidate.identifier() == StateVariable::Identifier::Custom) descriptors.push_back(candidate);

            // if there are custom state variables
            if (!descriptors.empty())
            {
                // create a text file
                TextOutFile out(this, itemName() + "_customstate_" + std::to_string(h), "custom state variables");

                // write the header
                out.addColumn("spatial cell index", "", 'd');
                for (const auto& descriptor : descriptors)
                {
                    string unit = descriptor.quantity();
                    if (!unit.empty()) unit = units->unit(unit);
                    out.addColumn(descriptor.description(), unit, descriptor.format());
                }

                // write a line for each cell
                for (int m = 0; m != numCells; ++m)
                {
                    vector<double> row({static_cast<double>(m)});
                    for (const auto& descriptor : descriptors)
                    {
                        double value = ms->custom(m, h, descriptor.customIndex());
                        if (!descriptor.quantity().empty()) value = units->out(descriptor.quantity(), value);
                        row.push_back(value);
                    }
                    out.writeRow(row);
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////
