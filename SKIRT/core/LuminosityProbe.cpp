/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LuminosityProbe.hpp"
#include "Configuration.hpp"
#include "SourceSystem.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void LuminosityProbe::probeSetup()
{
    auto units = find<Units>();

    // select "local" or default wavelength grid
    auto probeWavelengthGrid = find<Configuration>()->wavelengthGrid(wavelengthGrid());
    int numWavelengths = probeWavelengthGrid->numBins();

    // get the sources
    const auto& sources = find<SourceSystem>()->sources();
    int numSources = sources.size();

    // create a text file and add the columns
    TextOutFile file(this, probeName()+"_luminosities", "primary source luminosities");
    file.addColumn("wavelength", units->uwavelength());
    file.addColumn("specific luminosity", units->umonluminosity());
    file.addColumn("luminosity in bin", units->ubolluminosity());
    for (int i=0; i!=numSources; ++i)
        file.addColumn("luminosity fraction for source " + std::to_string(i+1));

    // write the rows
    for (int ell=0; ell!=numWavelengths; ++ell)
    {
        double lambda = probeWavelengthGrid->wavelength(ell);
        double dlambda = probeWavelengthGrid->effectiveWidth(ell);
        Array Llambdav(numSources);
        for (int i=0; i!=numSources; ++i) Llambdav[i] = sources[i]->specificLuminosity(lambda);
        double Llambdatot = Llambdav.sum();
        double Ltot = Llambdatot * dlambda;

        std::vector<double> row({units->owavelength(lambda),
                                 units->omonluminosityWavelength(lambda, Llambdatot),
                                 units->obolluminosity(Ltot)});
        for (int i=0; i!=numSources; ++i) row.push_back(Llambdatot > 0 ? Llambdav[i]/Llambdatot : 0);
        file.writeRow(row);
    }
}

////////////////////////////////////////////////////////////////////
