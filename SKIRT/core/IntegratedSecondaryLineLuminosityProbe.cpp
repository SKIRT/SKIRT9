/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "IntegratedSecondaryLineLuminosityProbe.hpp"
#include "Array.hpp"
#include "Configuration.hpp"
#include "Indices.hpp"
#include "LockFree.hpp"
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
    if (find<Configuration>()->hasGasEmission())
    {
        // locate the medium system and units system
        auto ms = find<MediumSystem>();
        int numCells = ms->numCells();
        auto units = find<Units>();

        // loop over medium components with line emission
        for (int h : ms->gasMediumIndices())
        {
            auto mix = ms->mix(0, h);
            if (mix->hasLineEmission())
            {
                // get the line centers
                const Array centers = mix->lineEmissionCenters();
                int numLines = centers.size();

                // accumulate the per-line luminosities across spatial cells (distributed across processes)
                Array perLineLum(numLines);
                find<ParallelFactory>()->parallelDistributed()->call(
                    numCells, [ms, h, numLines, &perLineLum](size_t firstIndex, size_t numIndices) {
                        for (size_t m = firstIndex; m != firstIndex + numIndices; ++m)
                        {
                            Array spectrum = ms->lineEmissionSpectrum(m, h);
                            // use lock-free addition to avoid race conditions
                            for (int k = 0; k != numLines; ++k) LockFree::add(perLineLum[k], spectrum[k]);
                        }
                    });
                ProcessManager::sumToRoot(perLineLum, true);

                // write a text file for this medium component
                TextOutFile file(this, itemName() + "_integrated_line_luminosities_" + std::to_string(h),
                                 "spatially integrated per-line luminosities for medium component "
                                     + std::to_string(h));
                file.writeLine("# Spatially integrated per-line luminosities for medium component "
                               + std::to_string(h));
                file.addColumn("wavelength; " + units->swavelength(), units->uwavelength());
                file.addColumn("luminosity", units->ubolluminosity());
                for (int k : Indices(centers, units->rwavelength()))
                {
                    file.writeRow({units->owavelength(centers[k]), units->obolluminosity(perLineLum[k])});
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////
