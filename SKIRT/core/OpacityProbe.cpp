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
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
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

        using MatType = MaterialMix::MaterialType;

        // system-level opacity
        auto systemCell = [numWaves, &wave, &ms](int m) {
            Array values(numWaves);
            for (int ell = 0; ell < numWaves; ++ell) values[ell] = ms->opacityExt(wave[ell], m);
            return values;
        };

        // type-level opacity
        auto typeCell = [numWaves, &wave, &ms](int m, MatType type) {
            Array values(numWaves);
            for (int ell = 0; ell < numWaves; ++ell) values[ell] = ms->opacityExt(wave[ell], m, type);
            return values;
        };

        // component-level opacity
        auto componentCell = [numWaves, &wave, &ms](int m, int h) {
            Array values(numWaves);
            for (int ell = 0; ell < numWaves; ++ell) values[ell] = ms->opacityExt(wave[ell], m, h);
            return values;
        };

        // fragment-level opacity
        auto fragmentCell = [numWaves, &wave, &ms](int m, int h, int f) {
            Array values(numWaves);
            auto mix = ms->mix(0, h)->find<FragmentDustMixDecorator>(false);
            for (int ell = 0; ell < numWaves; ++ell)
            {
                values[ell] = ms->callWithMaterialState(
                    [mix, f, lambda = wave[ell]](const MaterialState* mst) {
                        return mix->populationOpacityExt(f, lambda, mst);
                    },
                    m, h);
            }
            return values;
        };

        // define the call-back function to add column definitions
        auto addColumnDefinitions = [wlg, units](TextOutFile& outfile) {
            for (int ell : Indices(wlg->numBins(), units->rwavelength()))
            {
                outfile.addColumn("opacity at " + units->swavelength() + " = "
                                      + StringUtils::toString(units->owavelength(wlg->wavelength(ell)), 'g') + " "
                                      + units->uwavelength(),
                                  units->uopacity());
            }
        };

        ProbeFormBridge bridge(this, form());
        switch (aggregation())
        {
            case Aggregation::System:
            {
                bridge.writeQuantity("k", "tau", "opacity", "dimensionless", "opacity", "optical depth", axis,
                                     units->uwavelength(), addColumnDefinitions,
                                     [systemCell](int m) { return systemCell(m); });
                break;
            }

            case Aggregation::Type:
            {
                using MatType = MaterialMix::MaterialType;

                if (ms->hasDust())
                {
                    MatType type = MatType::Dust;
                    bridge.writeQuantity("dust_k", "dust_tau", "opacity", "dimensionless", "opacity", "optical depth",
                                         axis, units->uwavelength(), addColumnDefinitions,
                                         [typeCell, type](int m) { return typeCell(m, type); });
                }
                if (ms->hasElectrons())
                {
                    MatType type = MatType::Electrons;
                    bridge.writeQuantity("elec_k", "elec_tau", "opacity", "dimensionless", "opacity", "optical depth",
                                         axis, units->uwavelength(), addColumnDefinitions,
                                         [typeCell, type](int m) { return typeCell(m, type); });
                }
                if (ms->hasGas())
                {
                    MatType type = MatType::Gas;
                    bridge.writeQuantity("gas_k", "gas_tau", "opacity", "dimensionless", "opacity", "optical depth",
                                         axis, units->uwavelength(), addColumnDefinitions,
                                         [typeCell, type](int m) { return typeCell(m, type); });
                }
                break;
            }

            case Aggregation::Component:
            {
                int numMedia = ms->numMedia();
                for (int h = 0; h != numMedia; ++h)
                {
                    bridge.writeQuantity(std::to_string(h) + "_k", std::to_string(h) + "_tau", "opacity",
                                         "dimensionless", "opacity", "optical depth", axis, units->uwavelength(),
                                         addColumnDefinitions,
                                         [componentCell, h](int m) { return componentCell(m, h); });
                }
                break;
            }

            case Aggregation::Fragment:
            {
                for (int h : ms->dustMediumIndices())
                {
                    auto mix = ms->mix(0, h)->find<FragmentDustMixDecorator>(false);
                    if (!mix) continue;

                    int numPops = mix->numPopulations();
                    for (int f = 0; f != numPops; ++f)
                    {
                        bridge.writeQuantity(std::to_string(h) + "_" + std::to_string(f) + "_k",
                                             std::to_string(h) + "_" + std::to_string(f) + "_tau", "opacity",
                                             "dimensionless", "opacity", "optical depth", axis, units->uwavelength(),
                                             addColumnDefinitions,
                                             [fragmentCell, h, f](int m) { return fragmentCell(m, h, f); });
                    }
                }
                break;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////
