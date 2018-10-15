/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustTemperaturePerCellProbe.hpp"
#include "Configuration.hpp"
#include "MediumSystem.hpp"
#include "SpatialGrid.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void DustTemperaturePerCellProbe::probeRun()
{
    if (find<Configuration>()->hasRadiationField() && find<MediumSystem>()->hasDust())
    {
        // locate the medium system
        auto ms = find<MediumSystem>();
        auto grid = ms->grid();
        auto units = find<Units>();

        // create a text file
        TextOutFile file(this, itemName() + "_T", "dust temperature per cell");

        // write the header
        file.writeLine("# Indicative dust temperature per spatial cell");
        file.addColumn("spatial cell index", "", 'd');
        file.addColumn("indicative dust temperature", units->utemperature(), 'g');

        // write a line for each cell
        int numCells = grid->numCells();
        for (int m=0; m!=numCells; ++m)
        {
            const Array& Jv = ms->meanIntensity(m);
            double sumRhoT = 0.;
            double sumRho = 0.;
            for (int h=0; h!=ms->numMedia(); ++h)
            {
                if (ms->isDust(h))
                {
                    double rho = ms->massDensity(m,h);
                    if (rho > 0.)
                    {
                        double T = ms->mix(m,h)->equilibriumTemperature(Jv);
                        sumRhoT += rho*T;
                        sumRho += rho;
                    }
                }
            }
            double T = sumRho > 0. ? sumRhoT / sumRho : 0.;
            file.writeRow(m, units->otemperature(T));
        }
    }
}

////////////////////////////////////////////////////////////////////
