/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CustomStateProbe.hpp"
#include "Configuration.hpp"
#include "MediumSystem.hpp"
#include "NR.hpp"
#include "ProbeFormBridge.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // returns a list of integers corresponding to the indices and index ranges in the specified string
    // - if the string is empty, all indices 0...numIndices are included
    // - otherwise, the string must be a comma-separated list of zero-based indices or index ranges
    // portions of the string that don't conform to the syntax are ignored
    vector<int> listIndices(int n, string list)
    {
        vector<int> result;
        list = StringUtils::squeeze(list);
        if (list.empty())
        {
            for (int i = 0; i != n; ++i) result.push_back(i);
        }
        else
        {
            for (const string& segment : StringUtils::split(list, ","))
            {
                auto parts = StringUtils::split(segment, "-");
                if (parts.size() == 1)
                {
                    if (!StringUtils::isValidInt(parts[0])) continue;
                    int index = StringUtils::toInt(parts[0]);
                    if (index >= n) continue;
                    result.push_back(index);
                }
                else if (parts.size() == 2)
                {
                    if (!StringUtils::isValidInt(parts[0])) continue;
                    if (!StringUtils::isValidInt(parts[1])) continue;
                    int index1 = StringUtils::toInt(parts[0]);
                    int index2 = StringUtils::toInt(parts[1]);
                    if (index2 >= n) index2 = n - 1;
                    for (int i = index1; i <= index2; ++i) result.push_back(i);
                }
            }
        }
        return result;
    }
}

////////////////////////////////////////////////////////////////////

void CustomStateProbe::probe()
{
    if (find<Configuration>()->hasMedium())
    {
        // locate the medium system and units system
        auto ms = find<MediumSystem>();
        auto units = find<Units>();

        // construct a bridge
        ProbeFormBridge bridge(this, form());

        // loop over the media
        int numMedia = ms->numMedia();
        for (int h = 0; h != numMedia; ++h)
        {
            // get the custom state variable descriptors for this medium
            auto mix = ms->media()[h]->mix();
            vector<StateVariable> descriptors;
            for (const auto& candidate : mix->specificStateVariableInfo())
                if (candidate.identifier() == StateVariable::Identifier::Custom) descriptors.push_back(candidate);

            // get the indices corresponding to the user property and to the nr of custom variables for this medium
            auto indexv = listIndices(descriptors.size(), indices());

            // if there are custom variables to be probed for this medium
            if (!indexv.empty())
            {
                // construct the "quantity axis" from the zero-based indices
                Array axis = NR::array(indexv);

                // construct a list of corresponding unit conversion factors (for performance reasons)
                Array conv(axis.size());
                int outindex = 0;
                for (int index : indexv)
                {
                    const auto& descriptor = descriptors[index];
                    conv[outindex] = descriptor.quantity().empty() ? 1. : units->out(descriptor.quantity(), 1.);
                    outindex++;
                }

                // define the call-back function to add column definitions
                auto addColumnDefinitions = [&indexv, &descriptors, units](TextOutFile& outfile) {
                    for (int index : indexv)
                    {
                        const auto& descriptor = descriptors[index];
                        string qty = descriptor.quantity();
                        string unit = qty.empty() ? "1" : units->unit(qty);
                        outfile.addColumn(descriptor.description(), unit, descriptor.format());
                    }
                };

                // define the call-back function to retrieve a compound value in output units
                auto valueInCell = [ms, &indexv, &descriptors, h, conv](int m) {
                    Array out(indexv.size());
                    int outindex = 0;
                    for (int index : indexv)
                    {
                        const auto& descriptor = descriptors[index];
                        out[outindex] = conv[outindex] * ms->custom(m, h, descriptor.customIndex());
                        outindex++;
                    }
                    return out;
                };

                // define the call-back function to retrieve a cell weight; the units don't matter
                auto weightInCell = [ms](int m) { return ms->massDensity(m); };

                // get the units of the first variable being probed
                string qty = descriptors[indexv[0]].quantity();
                string unit = qty.empty() ? "1" : units->unit(qty);

                // output the requested variables for this medium
                string sh = std::to_string(h);
                bridge.writeQuantity(sh + "_customstate", unit, "custom state quantities",
                                     "density-weighted custom state quantities", axis, "1", addColumnDefinitions,
                                     valueInCell, weightInCell);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////
