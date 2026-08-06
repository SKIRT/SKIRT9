/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LuminosityProbe.hpp"
#include "Configuration.hpp"
#include "Indices.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "ProcessManager.hpp"
#include "SourceSystem.hpp"
#include "Table.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void LuminosityProbe::probe()
{
    auto units = find<Units>();

    // select "local" or default wavelength grid
    auto probeWavelengthGrid = find<Configuration>()->wavelengthGrid(wavelengthGrid());
    int numWavelengths = probeWavelengthGrid->numBins();

    // get the sources
    const auto& sources = find<SourceSystem>()->sources();
    int numSources = sources.size();

    // calculate the luminosities in parallel because this can take a while for imported sources with many entities
    Table<2> Llambdavv(numWavelengths, numSources);
    find<ParallelFactory>()->parallelDistributed()->call(
        numWavelengths, [&Llambdavv, probeWavelengthGrid, sources](size_t firstIndex, size_t numIndices) {
            int numSources = sources.size();
            for (size_t ell = firstIndex; ell != firstIndex + numIndices; ++ell)
            {
                double lambda = probeWavelengthGrid->wavelength(ell);
                for (int i = 0; i != numSources; ++i) Llambdavv(ell, i) = sources[i]->specificLuminosity(lambda);
            }
        });
    ProcessManager::sumToAll(Llambdavv.data());

    // create a text file and add the columns
    TextOutFile file(this, itemName() + "_luminosities", "primary source luminosities");
    file.addColumn("wavelength; " + units->swavelength(), units->uwavelength());
    file.addColumn("specific luminosity; " + units->smonluminosity(), units->umonluminosity());
    file.addColumn("luminosity in bin", units->ubolluminosity());
    for (int i = 0; i != numSources; ++i) file.addColumn("luminosity fraction for source " + std::to_string(i + 1));

    // write the rows
    for (int ell : Indices(numWavelengths, units->rwavelength()))
    {
        double lambda = probeWavelengthGrid->wavelength(ell);
        double dlambda = probeWavelengthGrid->effectiveWidth(ell);
        double Llambdatot = 0.;
        for (int i = 0; i != numSources; ++i) Llambdatot += Llambdavv(ell, i);
        double Ltot = Llambdatot * dlambda;

        std::vector<double> row(
            {units->owavelength(lambda), units->omonluminosity(lambda, Llambdatot), units->obolluminosity(Ltot)});
        for (int i = 0; i != numSources; ++i) row.push_back(Llambdatot > 0 ? Llambdavv(ell, i) / Llambdatot : 0);
        file.writeRow(row);
    }
}

////////////////////////////////////////////////////////////////////
