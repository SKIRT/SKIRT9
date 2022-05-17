/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MetallicityProbe.hpp"
#include "Configuration.hpp"
#include "MediumSystem.hpp"
#include "ProbeFormBridge.hpp"

////////////////////////////////////////////////////////////////////

void MetallicityProbe::probe()
{
    if (find<Configuration>()->hasMedium())
    {
        // locate the medium system
        auto ms = find<MediumSystem>();

        // build a list of indices for gas medium components that store metallicity
        vector<int> hv;
        for (int h : ms->gasMediumIndices())
        {
            auto mix = ms->media()[h]->mix();
            for (const auto& candidate : mix->specificStateVariableInfo())
            {
                if (candidate.identifier() == StateVariable::Identifier::Metallicity)
                {
                    hv.push_back(h);
                    break;
                }
            }
        }

        // construct a bridge
        ProbeFormBridge bridge(this, form());

        // produce output for each of the detected components, if any
        for (int h : hv)
        {
            string sh = std::to_string(h);

            // define the call-back functions that will be passed to the bridge
            auto valueInCell = [ms, h](int m) { return ms->metallicity(m, h); };
            auto weightInCell = [ms, h](int m) { return ms->massDensity(m, h); };

            // produce output
            bridge.writeQuantity(sh + "_Z", "dimensionless", "metallicity", "density-weighted metallicity", valueInCell,
                                 weightInCell);
        }
    }
}

////////////////////////////////////////////////////////////////////
