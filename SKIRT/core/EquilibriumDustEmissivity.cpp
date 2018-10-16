/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "EquilibriumDustEmissivity.hpp"
#include "Configuration.hpp"
#include "MaterialMix.hpp"
#include "PlanckFunction.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

Array EquilibriumDustEmissivity::emissivity(const MaterialMix* mix, const Array& Jv) const
{
    // get the output wavelength grid
    auto wavelengthGrid = find<Configuration>()->emissionSpectrumWLG();
    int numWavelengths = wavelengthGrid->numBins();

    // accumulate the emissivities at the equilibrium temperature for all populations in the dust mix
    Array ev(numWavelengths);
    {
        double T = mix->equilibriumTemperature(Jv);
        PlanckFunction B(T);
        for (int ell=0; ell<numWavelengths; ell++)
        {
            double lambda = wavelengthGrid-> wavelength(ell);
            ev[ell] += mix->sectionAbs(lambda) * B(lambda);
        }
    }

    // convert emissivity from "per hydrogen atom" to "per unit mass"
    ev /= mix->mass();
    return ev;
}

////////////////////////////////////////////////////////////////////
