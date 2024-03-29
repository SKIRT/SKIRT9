/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DensityProbe.hpp"
#include "Configuration.hpp"
#include "FragmentDustMixDecorator.hpp"
#include "MediumSystem.hpp"
#include "ProbeFormBridge.hpp"

////////////////////////////////////////////////////////////////////

void DensityProbe::probe()
{
    if (find<Configuration>()->hasMedium())
    {
        // locate the medium system
        auto ms = find<MediumSystem>();

        // construct a bridge
        ProbeFormBridge bridge(this, form());

        // produce output depending on aggregation level
        switch (aggregation())
        {
            case Aggregation::Type:
            {
                if (ms->hasDust())
                {
                    bridge.writeQuantity("dust_rho", "dust_Sigma", "massvolumedensity", "masssurfacedensity",
                                         "mass density", "mass surface density",
                                         [ms](int m) { return ms->dustMassDensity(m); });
                }
                if (ms->hasElectrons())
                {
                    bridge.writeQuantity("elec_n", "elec_N", "numbervolumedensity", "numbersurfacedensity",
                                         "number density", "column density",
                                         [ms](int m) { return ms->electronNumberDensity(m); });
                }
                if (ms->hasGas())
                {
                    bridge.writeQuantity("gas_n", "gas_N", "numbervolumedensity", "numbersurfacedensity",
                                         "number density", "column density",
                                         [ms](int m) { return ms->gasNumberDensity(m); });
                }
                break;
            }
            case Aggregation::Component:
            {
                int numMedia = ms->numMedia();
                for (int h = 0; h != numMedia; ++h)
                {
                    string sh = std::to_string(h);
                    if (ms->isDust(h))
                    {
                        bridge.writeQuantity(sh + "_rho", sh + "_Sigma", "massvolumedensity", "masssurfacedensity",
                                             "mass density", "mass surface density",
                                             [ms, h](int m) { return ms->massDensity(m, h); });
                    }
                    else
                    {
                        bridge.writeQuantity(sh + "_n", sh + "_N", "numbervolumedensity", "numbersurfacedensity",
                                             "number density", "column density",
                                             [ms, h](int m) { return ms->numberDensity(m, h); });
                    }
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
                            auto densityInCell = [ms, mix, h, f](int m) {
                                return mix->populationMass(f)
                                       * ms->callWithMaterialState(
                                           [mix, f](const MaterialState* mst) {
                                               return mix->populationNumberDensity(f, mst);
                                           },
                                           m, h);
                            };

                            string shf = std::to_string(h) + "_" + std::to_string(f);
                            bridge.writeQuantity(shf + "_rho", shf + "_Sigma", "massvolumedensity",
                                                 "masssurfacedensity", "mass density", "mass surface density",
                                                 densityInCell);
                        }
                    }
                }
                break;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////
