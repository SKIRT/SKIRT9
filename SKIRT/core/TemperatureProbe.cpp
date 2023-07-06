/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TemperatureProbe.hpp"
#include "Configuration.hpp"
#include "FragmentDustMixDecorator.hpp"
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
        // can probe dust temperature only during or after secondary emission phase
        // issue a warning because otherwise this is a very confusing problem for the user to locate and fix
        if (probeAfter() == ProbeAfter::Setup || probeAfter() == ProbeAfter::Primary)
        {
            find<Log>()->warning(
                typeAndName()
                + " can probe dust temperature only during or after secondary emission; adjust probeAfter property");
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
                case Aggregation::Fragment:
                {
                    for (int h : ms->dustMediumIndices())
                    {
                        auto mix = ms->mix(0, h)->find<FragmentDustMixDecorator>(false);
                        if (mix)
                        {
                            int numPops = mix->numPopulations();
                            for (int f = 0; f != numPops; ++f)
                            {
                                auto weightInCell = [ms, mix, h, f](int m) {
                                    return ms->callWithMaterialState(
                                        [mix, f](const MaterialState* mst) {
                                            return mix->populationNumberDensity(f, mst);
                                        },
                                        m, h);
                                };
                                auto valueInCell = [ms, mix, f, &weightInCell](int m) {
                                    if (weightInCell(m) <= 0.) return 0.;
                                    const Array& Jv = ms->meanIntensity(m);
                                    return mix->populationTemperature(f, Jv);
                                };

                                string shf = std::to_string(h) + "_" + std::to_string(f);
                                bridge.writeQuantity(shf + "_T", "temperature", "indicative temperature",
                                                     "density-weighted indicative temperature", valueInCell,
                                                     weightInCell);
                            }
                        }
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
            case Aggregation::Fragment: break;
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
            case Aggregation::Fragment: break;
        }
    }
}

////////////////////////////////////////////////////////////////////
