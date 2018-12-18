/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "EquilibriumDustEmissivity.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "MultiGrainDustMix.hpp"
#include "PlanckFunction.hpp"

////////////////////////////////////////////////////////////////////

Array EquilibriumDustEmissivity::emissivity(const MaterialMix* mix, const Array& Jv) const
{
    // get the output wavelength grid
    auto wavelengthGrid = find<Configuration>()->dustEmissionWLG();
    int numWavelengths = wavelengthGrid->numBins();

    // determine the emissivity spectrum
    Array ev(numWavelengths);
    {
        // if the mix exposes multiple representative dust grains (size bins),
        // accumulate the emissivities at the equilibrium temperature for all representative grains
        auto mgmix = useSingleGrain() ? nullptr : dynamic_cast<const MultiGrainDustMix*>(mix);
        if (mgmix)
        {
            int numBins = mgmix->numBins();
            for (int b=0; b!=numBins; ++b)
            {
                double T = mgmix->binEquilibriumTemperature(b, Jv);
                PlanckFunction B(T);
                for (int ell=0; ell<numWavelengths; ell++)
                {
                    double lambda = wavelengthGrid-> wavelength(ell);
                    ev[ell] += mgmix->binSectionAbs(b, lambda) * B(lambda);
                }
            }
        }

        // otherwise, just use a single representative grain for the complete mix
        else
        {
            double T = mix->equilibriumTemperature(Jv);
            PlanckFunction B(T);
            for (int ell=0; ell<numWavelengths; ell++)
            {
                double lambda = wavelengthGrid-> wavelength(ell);
                ev[ell] += mix->sectionAbs(lambda) * B(lambda);
            }
        }
    }

    // convert emissivity from "per hydrogen atom" to "per unit mass"
    ev /= mix->mass();
    return ev;
}

////////////////////////////////////////////////////////////////////
