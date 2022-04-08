/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TemperatureProbe.hpp"
#include "Configuration.hpp"
#include "MediumSystem.hpp"
#include "ProbeFormBridge.hpp"

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

    // only if probing after the full simulation run and if panchromatic radiation field is available
    if (probeAfter() == ProbeAfter::Run && config->hasPanRadiationField() && ms->hasDust())
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

    // ------- handle electrons -------

    if (ms->hasElectrons())
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

    if (ms->hasGas())
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
