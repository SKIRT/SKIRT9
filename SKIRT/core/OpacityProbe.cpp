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
        // locate the medium system and construct lists of wavelength information in **output order**
        auto units = find<Units>();
        auto ms = find<MediumSystem>();
        auto wlg = find<Configuration>()->wavelengthGrid(wavelengthGrid());
        int numWaves = wlg->numBins();
        Array wave(numWaves);  // wavelength in internal units
        Array axis(numWaves);  // wavelength in output units
        int outell = 0;
        for (int ell : Indices(numWaves, units->rwavelength()))
        {
            wave[outell] = wlg->wavelength(ell);
            axis[outell] = units->owavelength(wlg->wavelength(ell));
            outell++;
        }

        // define the call-back function to add column definitions
        auto addColumnDefinitions = [axis, units](TextOutFile& outfile) {
            for (double outwave : axis)
            {
                outfile.addColumn("opacity at " + units->swavelength() + " = " + StringUtils::toString(outwave, 'g')
                                      + " " + units->uwavelength(),
                                  units->uopacity());
            }
        };

        // construct a bridge and tell it to produce output of the user-configured kind
        ProbeFormBridge bridge(this, form());
        switch (aggregation())
        {
            case Aggregation::System:
            {
                auto valueInCell = [numWaves, &wave, &ms](int m) {
                    Array values(numWaves);
                    for (int ell = 0; ell < numWaves; ++ell) values[ell] = ms->opacityExt(wave[ell], m);
                    return values;
                };

                bridge.writeQuantity("k", "tau", "opacity", "dimensionless", "opacity", "optical depth", axis,
                                     units->uwavelength(), addColumnDefinitions, valueInCell);
                break;
            }

            case Aggregation::Type:
            {
                using MatType = MaterialMix::MaterialType;

                auto valueInCell = [numWaves, &wave, &ms](int m, MatType type) {
                    Array values(numWaves);
                    for (int ell = 0; ell < numWaves; ++ell) values[ell] = ms->opacityExt(wave[ell], m, type);
                    return values;
                };

                if (ms->hasDust())
                {
                    bridge.writeQuantity("dust_k", "dust_tau", "opacity", "dimensionless", "opacity", "optical depth",
                                         axis, units->uwavelength(), addColumnDefinitions,
                                         [valueInCell](int m) { return valueInCell(m, MatType::Dust); });
                }
                if (ms->hasElectrons())
                {
                    bridge.writeQuantity("elec_k", "elec_tau", "opacity", "dimensionless", "opacity", "optical depth",
                                         axis, units->uwavelength(), addColumnDefinitions,
                                         [valueInCell](int m) { return valueInCell(m, MatType::Electrons); });
                }
                if (ms->hasGas())
                {
                    bridge.writeQuantity("gas_k", "gas_tau", "opacity", "dimensionless", "opacity", "optical depth",
                                         axis, units->uwavelength(), addColumnDefinitions,
                                         [valueInCell](int m) { return valueInCell(m, MatType::Gas); });
                }
                break;
            }

            case Aggregation::Component:
            {
                auto valueInCell = [numWaves, &wave, &ms](int m, int h) {
                    Array values(numWaves);
                    for (int ell = 0; ell < numWaves; ++ell) values[ell] = ms->opacityExt(wave[ell], m, h);
                    return values;
                };

                int numMedia = ms->numMedia();
                for (int h = 0; h != numMedia; ++h)
                {
                    bridge.writeQuantity(std::to_string(h) + "_k", std::to_string(h) + "_tau", "opacity",
                                         "dimensionless", "opacity", "optical depth", axis, units->uwavelength(),
                                         addColumnDefinitions, [valueInCell, h](int m) { return valueInCell(m, h); });
                }
                break;
            }

            case Aggregation::Fragment:
            {
                auto valueInCell = [numWaves, &wave, &ms](int m, int h, int f) {
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

                for (int h : ms->dustMediumIndices())
                {
                    auto mix = ms->mix(0, h)->find<FragmentDustMixDecorator>(false);
                    if (mix)
                    {
                        int numPops = mix->numPopulations();
                        for (int f = 0; f != numPops; ++f)
                        {
                            bridge.writeQuantity(std::to_string(h) + "_" + std::to_string(f) + "_k",
                                                 std::to_string(h) + "_" + std::to_string(f) + "_tau", "opacity",
                                                 "dimensionless", "opacity", "optical depth", axis,
                                                 units->uwavelength(), addColumnDefinitions,
                                                 [valueInCell, h, f](int m) { return valueInCell(m, h, f); });
                        }
                    }
                }
                break;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////
