/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ConvergenceCutsProbe.hpp"
#include "Configuration.hpp"
#include "DefaultCutsForm.hpp"
#include "MediumSystem.hpp"
#include "ProbeFormBridge.hpp"

////////////////////////////////////////////////////////////////////

void ConvergenceCutsProbe::setupSelfBefore()
{
    SpecialtyWhenProbe::setupSelfBefore();

    _form = new DefaultCutsForm(this);
}

////////////////////////////////////////////////////////////////////

void ConvergenceCutsProbe::probe()
{
    if (find<Configuration>()->hasMedium())
    {
        // locate the medium system
        auto ms = find<MediumSystem>();

        // construct a bridge with our hard-coded form
        ProbeFormBridge bridge(this, _form);

        // write true and gridded density cuts
        // (we don't provide all info for projected quantities because we know our form will never use these)
        if (ms->hasDust())
        {
            bridge.writeQuantity(
                "dust_t", "", "massvolumedensity", "masssurfacedensity", "mass density", "",
                [ms](Position bfr) { return ms->dustMassDensity(bfr); }, nullptr);
            bridge.writeQuantity("dust_g", "", "massvolumedensity", "masssurfacedensity", "mass density", "",
                                 [ms](int m) { return ms->dustMassDensity(m); });
        }
        if (ms->hasElectrons())
        {
            bridge.writeQuantity(
                "elec_t", "", "numbervolumedensity", "numbersurfacedensity", "number density", "",
                [ms](Position bfr) { return ms->electronNumberDensity(bfr); }, nullptr);
            bridge.writeQuantity("elec_g", "", "numbervolumedensity", "numbersurfacedensity", "number density", "",
                                 [ms](int m) { return ms->electronNumberDensity(m); });
        }
        if (ms->hasGas())
        {
            bridge.writeQuantity(
                "gas_t", "", "numbervolumedensity", "numbersurfacedensity", "number density", "",
                [ms](Position bfr) { return ms->gasNumberDensity(bfr); }, nullptr);
            bridge.writeQuantity("gas_g", "", "numbervolumedensity", "numbersurfacedensity", "number density", "",
                                 [ms](int m) { return ms->gasNumberDensity(m); });
        }
    }
}

////////////////////////////////////////////////////////////////////
