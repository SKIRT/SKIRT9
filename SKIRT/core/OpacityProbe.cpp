/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "OpacityProbe.hpp"
#include "Configuration.hpp"
#include "FragmentDustMixDecorator.hpp"
#include "MediumSystem.hpp"
#include "ProbeFormBridge.hpp"

////////////////////////////////////////////////////////////////////

Range OpacityProbe::wavelengthRange() const
{
    return Range(wavelength(), wavelength());
}

////////////////////////////////////////////////////////////////////

void OpacityProbe::probe()
{
    if (find<Configuration>()->hasMedium())
    {
        // locate the medium system and cache the wavelength
        auto ms = find<MediumSystem>();
        double lambda = wavelength();

        // construct a bridge
        ProbeFormBridge bridge(this, form());

        // produce output depending on aggregation level
        switch (aggregation())
        {
            case Aggregation::System:
            {
                bridge.writeQuantity("k", "tau", "opacity", "dimensionless", "opacity", "optical depth",
                                     [ms, lambda](int m) { return ms->opacityExt(lambda, m); });
                break;
            }
            case Aggregation::Type:
            {
                using MatType = MaterialMix::MaterialType;
                if (ms->hasDust())
                {
                    bridge.writeQuantity("dust_k", "dust_tau", "opacity", "dimensionless", "opacity", "optical depth",
                                         [ms, lambda](int m) { return ms->opacityExt(lambda, m, MatType::Dust); });
                }
                if (ms->hasElectrons())
                {
                    bridge.writeQuantity("elec_k", "elec_tau", "opacity", "dimensionless", "opacity", "optical depth",
                                         [ms, lambda](int m) { return ms->opacityExt(lambda, m, MatType::Electrons); });
                }
                if (ms->hasGas())
                {
                    bridge.writeQuantity("gas_k", "gas_tau", "opacity", "dimensionless", "opacity", "optical depth",
                                         [ms, lambda](int m) { return ms->opacityExt(lambda, m, MatType::Gas); });
                }
                break;
            }
            case Aggregation::Component:
            {
                int numMedia = ms->numMedia();
                for (int h = 0; h != numMedia; ++h)
                {
                    string sh = std::to_string(h);
                    bridge.writeQuantity(sh + "_k", sh + "_tau", "opacity", "dimensionless", "opacity", "optical depth",
                                         [ms, lambda, h](int m) { return ms->opacityExt(lambda, m, h); });
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
                            auto opacityInCell = [ms, lambda, mix, h, f](int m) {
                                return ms->callWithMaterialState(
                                    [lambda, mix, f](const MaterialState* mst) {
                                        return mix->populationOpacityExt(f, lambda, mst);
                                    },
                                    m, h);
                            };

                            string shf = std::to_string(h) + "_" + std::to_string(f);
                            bridge.writeQuantity(shf + "_k", shf + "_tau", "opacity", "dimensionless", "opacity",
                                                 "optical depth", opacityInCell);
                        }
                    }
                }
                break;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////
