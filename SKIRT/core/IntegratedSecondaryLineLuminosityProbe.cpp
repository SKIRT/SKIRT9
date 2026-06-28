/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "IntegratedSecondaryLineLuminosityProbe.hpp"
#include "Array.hpp"
#include "MaterialMix.hpp"
#include "MediumSystem.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "ProcessManager.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void IntegratedSecondaryLineLuminosityProbe::probe()
{
    auto ms = find<MediumSystem>();
    if (!ms) return;

    auto units = find<Units>();
    const int numMedia = ms->numMedia();
    const int numCells = ms->numCells();

    // visit each medium component and accumulate per-line totals when its mix supports line emission
    for (int h = 0; h != numMedia; ++h)
    {
        const auto* mix = ms->mix(0, h);
        if (!mix || !mix->hasLineEmission()) continue;

        const Array centers = mix->lineEmissionCenters();
        const int numLines = centers.size();
        if (numLines == 0) continue;

        // accumulate the per-line luminosities across spatial cells (distributed across processes)
        Array perLineLum(numLines);
        find<ParallelFactory>()->parallelDistributed()->call(
            numCells, [ms, h, &perLineLum](size_t firstIndex, size_t numIndices) {
                for (size_t m = firstIndex; m != firstIndex + numIndices; ++m)
                {
                    Array spectrum = ms->lineEmissionSpectrum(m, h);
                    if (spectrum.size() == perLineLum.size()) perLineLum += spectrum;
                }
            });
        ProcessManager::sumToAll(perLineLum);

        // write a text file for this component
        TextOutFile file(this, itemName() + "_integrated_line_luminosities_" + std::to_string(h),
                         "spatially integrated per-line luminosities for medium component " + std::to_string(h));
        file.addColumn("wavelength; " + units->swavelength(), units->uwavelength());
        file.addColumn("luminosity", units->ubolluminosity());
        for (int k = 0; k != numLines; ++k)
        {
            file.writeRow({units->owavelength(centers[k]), units->obolluminosity(perLineLum[k])});
        }
    }
}

////////////////////////////////////////////////////////////////////
