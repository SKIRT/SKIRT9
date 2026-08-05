/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpatialCellPropertiesProbe.hpp"
#include "MediumSystem.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void SpatialCellPropertiesProbe::probe()
{
    // locate the medium system and the grid (it is OK for the medium system to have no media components)
    auto ms = find<MediumSystem>(false);
    if (ms)
    {
        auto grid = ms->grid();
        auto units = find<Units>();

        // create a text file
        TextOutFile out(this, itemName() + "_cellprops", "spatial cell properties");

        // write the header
        out.addColumn("spatial cell index", "", 'd');
        out.addColumn("x coordinate of cell center", units->ulength());
        out.addColumn("y coordinate of cell center", units->ulength());
        out.addColumn("z coordinate of cell center", units->ulength());
        out.addColumn("cell volume", units->uvolume());
        out.addColumn("optical depth of cell diagonal at "
                      + StringUtils::toString(units->owavelength(wavelength()), 'g') + " " + units->uwavelength());
        out.addColumn("dust mass density in cell", units->umassvolumedensity());
        out.addColumn("electron number density in cell", units->unumbervolumedensity());
        out.addColumn("hydrogen number density in cell", units->unumbervolumedensity());

        // write a line for each cell
        int numMedia = ms->numMedia();
        int numCells = grid->numCells();
        for (int m = 0; m != numCells; ++m)
        {
            Position p = grid->centralPositionInCell(m);
            double V = ms->volume(m);
            double tau = grid->diagonal(m) * ms->opacityExt(wavelength(), m);
            double dust = 0.;
            double elec = 0.;
            double gas = 0.;
            for (int h = 0; h != numMedia; ++h)
            {
                if (ms->isDust(h)) dust += ms->massDensity(m, h);
                if (ms->isElectrons(h)) elec += ms->numberDensity(m, h);
                if (ms->isGas(h)) gas += ms->numberDensity(m, h);
            }
            out.writeRow(vector<double>({static_cast<double>(m), units->olength(p.x()), units->olength(p.y()),
                                         units->olength(p.z()), units->ovolume(V), tau, units->omassvolumedensity(dust),
                                         units->onumbervolumedensity(elec), units->onumbervolumedensity(gas)}));
        }
    }
}

////////////////////////////////////////////////////////////////////
