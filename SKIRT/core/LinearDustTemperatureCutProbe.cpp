/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LinearDustTemperatureCutProbe.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "MediumSystem.hpp"
#include "SpatialGrid.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void LinearDustTemperatureCutProbe::probeRun()
{
    if (find<Configuration>()->hasRadiationField() && find<MediumSystem>()->hasDust())
    {
        // locate the medium system
        auto ms = find<MediumSystem>();
        auto grid = ms->grid();
        auto units = find<Units>();

        // get a characteristic size of the spatial grid
        double size = grid->boundingBox().diagonal();

        // get the line segment parameters
        Position p1(_startX, _startY, _startZ);
        Position p2(_endX, _endY, _endZ);
        double length = (p2-p1).norm();
        if (length < size*1e-10) throw FATALERROR("Line segment is too short");

        // create a text file
        TextOutFile file(this, itemName() + "_T", "dust temperature along a line segment");

        // write the header
        file.writeLine("# Indicative dust temperature along a given line segment");
        file.addColumn("distance from starting point", units->ulength());
        file.addColumn("indicative dust temperature", units->utemperature(), 'g');

        // write a line for each sample
        for (int i=0; i!=_numSamples; ++i)
        {
            // determine the sample position and the distance along the line
            double fraction = static_cast<double>(i) / static_cast<double>(_numSamples-1);
            Position p(p1 + fraction*(p2-p1));
            double distance = (p-p1).norm();

            // determine the cell containing the sample position
            int m = grid->cellIndex(p);

            // calculate the indicative dust temperature in that cell
            double Tout = 0.;
            if (m>=0)
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
                if (sumRho > 0.) Tout = sumRhoT / sumRho;
            }

            // write the row
            file.writeRow(units->olength(distance), units->otemperature(Tout));
        }
    }
}

////////////////////////////////////////////////////////////////////
