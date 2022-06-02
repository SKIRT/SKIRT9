/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TemperatureProbe.hpp"
#include "Configuration.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "ProbeFormBridge.hpp"

////////////////////////////////////////////////////////////////////

// return true if any of the medium components with the given indices store temperature in the medium state
bool hasTemperature(const MediumSystem* ms, const vector<int>& hv)
{
    for (int h : hv)
    {
        auto mix = ms->media()[h]->mix();
        for (const auto& candidate : mix->specificStateVariableInfo())
        {
            if (candidate.identifier() == StateVariable::Identifier::Temperature) return true;
        }
    }
    return false;
}

////////////////////////////////////////////////////////////////////

void TemperatureProbe::probe()
{
    auto config = find<Configuration>();

    // locate the medium system, or exit if there is none
    if (!config->hasMedium()) return;
    auto ms = find<MediumSystem>();

    // construct a bridge
    ProbeFormBridge bridge(this, form());

    // ------- handle dust -------

    // can probe dust temperature only if panchromatic radiation field is available
    if (config->hasPanRadiationField() && ms->hasDust())
    {
        // can probe dust temperature only after the full simulation run
        if (probeAfter() == ProbeAfter::Setup)
        {
            // issue a warning because (1) the default value is Setup and (2) otherwise this is a very confusing problem
            find<Log>()->warning("Cannot probe dust temperature after setup; set probeAfter to 'Run'");
        }
        else
        {
            // produce output depending on aggregation level
            switch (aggregation())
            {
                case Aggregation::Type:
                {
                    auto valueInCell = [ms](int m) { return ms->indicativeDustTemperature(m); };
                    auto weightInCell = [ms](int m) { return ms->dustMassDensity(m); };
                    bridge.writeQuantity("dust_T", "temperature", "indicative temperature",
                                         "density-weighted indicative temperature", valueInCell, weightInCell);
                    break;
                }
                case Aggregation::Component:
                {
                    for (int h : ms->dustMediumIndices())
                    {
                        string sh = std::to_string(h);
                        auto valueInCell = [ms, h](int m) { return ms->indicativeTemperature(m, h); };
                        auto weightInCell = [ms, h](int m) { return ms->massDensity(m, h); };
                        bridge.writeQuantity(sh + "_T", "temperature", "indicative temperature",
                                             "density-weighted indicative temperature", valueInCell, weightInCell);
                    }
                    break;
                }
            }
        }
    }

    // ------- handle electrons -------

    if (hasTemperature(ms, ms->electronMediumIndices()))
    {
        // produce output depending on aggregation level
        switch (aggregation())
        {
            case Aggregation::Type:
            {
                auto valueInCell = [ms](int m) { return ms->indicativeElectronTemperature(m); };
                auto weightInCell = [ms](int m) { return ms->electronNumberDensity(m); };
                bridge.writeQuantity("elec_T", "temperature", "indicative temperature",
                                     "density-weighted indicative temperature", valueInCell, weightInCell);
                break;
            }
            case Aggregation::Component:
            {
                for (int h : ms->electronMediumIndices())
                {
                    string sh = std::to_string(h);
                    auto valueInCell = [ms, h](int m) { return ms->indicativeTemperature(m, h); };
                    auto weightInCell = [ms, h](int m) { return ms->numberDensity(m, h); };
                    bridge.writeQuantity(sh + "_T", "temperature", "indicative temperature",
                                         "density-weighted indicative temperature", valueInCell, weightInCell);
                }
                break;
            }
        }
    }

    // ------- handle gas -------

    if (hasTemperature(ms, ms->gasMediumIndices()))
    {
        // produce output depending on aggregation level
        switch (aggregation())
        {
            case Aggregation::Type:
            {
                auto valueInCell = [ms](int m) { return ms->indicativeGasTemperature(m); };
                auto weightInCell = [ms](int m) { return ms->gasNumberDensity(m); };
                bridge.writeQuantity("gas_T", "temperature", "indicative temperature",
                                     "density-weighted indicative temperature", valueInCell, weightInCell);
                break;
            }
            case Aggregation::Component:
            {
                for (int h : ms->gasMediumIndices())
                {
                    string sh = std::to_string(h);
                    auto valueInCell = [ms, h](int m) { return ms->indicativeTemperature(m, h); };
                    auto weightInCell = [ms, h](int m) { return ms->numberDensity(m, h); };
                    bridge.writeQuantity(sh + "_T", "temperature", "indicative temperature",
                                         "density-weighted indicative temperature", valueInCell, weightInCell);
                }
                break;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////
