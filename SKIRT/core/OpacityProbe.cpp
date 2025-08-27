/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "OpacityProbe.hpp"
#include "Configuration.hpp"
#include "FragmentDustMixDecorator.hpp"
#include "Indices.hpp"
#include "MediumSystem.hpp"
#include "ProbeFormBridge.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

Range OpacityProbe::wavelengthRange() const
{
    if (wavelengthGrid())
    {
        wavelengthGrid()->setup();
        return wavelengthGrid()->wavelengthRange();
    }
    return Range();
}

////////////////////////////////////////////////////////////////////

WavelengthGrid* OpacityProbe::materialWavelengthGrid() const
{
    return wavelengthGrid();
}

////////////////////////////////////////////////////////////////////

void OpacityProbe::probe()
{
    if (find<Configuration>()->hasMedium())
    {
        // locate the medium system and cache the wavelength
        auto units = find<Units>();
        auto ms = find<MediumSystem>();
        auto wlg = find<Configuration>()->wavelengthGrid(wavelengthGrid());
        int numWaves = wlg->numBins();
        Array wave(numWaves);  // wavelength in internal units
        Array axis(numWaves);  // wavelength in output units
        for (int ell : Indices(numWaves, units->rwavelength()))
        {
            wave[ell] = wlg->wavelength(ell);
            axis[ell] = units->owavelength(wlg->wavelength(ell));
        }

        // use std::optional? or default values?
        using MatType = MaterialMix::MaterialType;
        auto valueInCell = [numWaves, &wave, &ms](int m, Aggregation agg, MatType* type, int* h, int* f) {
            Array values(numWaves);
            for (int ell = 0; ell < numWaves; ++ell)
            {
                switch (agg)
                {
                    case Aggregation::System: values[ell] = ms->opacityExt(wave[ell], m); break;
                    case Aggregation::Type: values[ell] = ms->opacityExt(wave[ell], m, *type); break;
                    case Aggregation::Component: values[ell] = ms->opacityExt(wave[ell], m, *h); break;
                    case Aggregation::Fragment:
                        auto mix = ms->mix(0, *h)->find<FragmentDustMixDecorator>(false);
                        values[ell] = ms->callWithMaterialState(
                            [mix, f, lambda = wave[ell]](const MaterialState* mst) {
                                return mix->populationOpacityExt(*f, lambda, mst);
                            },
                            m, *h);
                        break;
                }
            }
            return values;
        };

        ProbeFormBridge bridge(this, form());
        // produce output depending on aggregation level
        switch (aggregation())
        {
            case Aggregation::System:
            {
                auto systemValue = [valueInCell](int m) {
                    return valueInCell(m, Aggregation::System, nullptr, nullptr, nullptr);
                };

                bridge.writeQuantity("k", "tau", "opacity", "dimensionless", "opacity", "optical depth", axis,
                                     units->uwavelength(), nullptr, systemValue);
                break;
            }
            case Aggregation::Type:
            {
                using MatType = MaterialMix::MaterialType;
                if (ms->hasDust())
                {
                    MatType type = MatType::Dust;
                    auto dustValue = [valueInCell, &type](int m) {
                        return valueInCell(m, Aggregation::Type, &type, nullptr, nullptr);
                    };

                    bridge.writeQuantity("dust_k", "dust_tau", "opacity", "dimensionless", "opacity", "optical depth",
                                         axis, units->uwavelength(), nullptr, dustValue);
                }
                if (ms->hasElectrons())
                {
                    MatType type = MatType::Electrons;
                    auto elecValue = [valueInCell, &type](int m) {
                        return valueInCell(m, Aggregation::Type, &type, nullptr, nullptr);
                    };

                    bridge.writeQuantity("elec_k", "elec_tau", "opacity", "dimensionless", "opacity", "optical depth",
                                         axis, units->uwavelength(), nullptr, elecValue);
                }
                if (ms->hasGas())
                {
                    MatType type = MatType::Gas;
                    auto gasValue = [valueInCell, &type](int m) {
                        return valueInCell(m, Aggregation::Type, &type, nullptr, nullptr);
                    };

                    bridge.writeQuantity("gas_k", "gas_tau", "opacity", "dimensionless", "opacity", "optical depth",
                                         axis, units->uwavelength(), nullptr, gasValue);
                }
                break;
            }
            case Aggregation::Component:
            {
                int numMedia = ms->numMedia();
                for (int h = 0; h != numMedia; ++h)
                {
                    auto compValue = [valueInCell, &h](int m) {
                        return valueInCell(m, Aggregation::Component, nullptr, &h, nullptr);
                    };

                    string sh = std::to_string(h);
                    bridge.writeQuantity(sh + "_k", sh + "_tau", "opacity", "dimensionless", "opacity", "optical depth",
                                         axis, units->uwavelength(), nullptr, compValue);
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
                            auto fragValue = [valueInCell, &h, &f](int m) {
                                return valueInCell(m, Aggregation::Fragment, nullptr, &h, &f);
                            };

                            string shf = std::to_string(h) + "_" + std::to_string(f);
                            bridge.writeQuantity(shf + "_k", shf + "_tau", "opacity", "dimensionless", "opacity",
                                                 "optical depth", axis, units->uwavelength(), nullptr, fragValue);
                        }
                    }
                }
                break;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////
