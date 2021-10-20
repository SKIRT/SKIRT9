/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PlanarGasTemperatureCutsProbe.hpp"
#include "Array.hpp"
#include "Configuration.hpp"
#include "FITSInOut.hpp"
#include "MediumSystem.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void PlanarGasTemperatureCutsProbe::writeGasTemperatureCut(Probe* probe, bool xd, bool yd, bool zd, double xc,
                                                           double yc, double zc, int Nx, int Ny, int Nz)
{
    // locate relevant simulation items
    auto units = probe->find<Units>();
    auto ms = probe->find<MediumSystem>();
    auto grid = ms->grid();

    // determine spatial configuration (regardless of cut direction)
    Box box = grid->boundingBox();
    double xpsize = box.xwidth() / Nx;
    double ypsize = box.ywidth() / Ny;
    double zpsize = box.zwidth() / Nz;
    double xbase = box.xmin() + 0.5 * xpsize;
    double ybase = box.ymin() + 0.5 * ypsize;
    double zbase = box.zmin() + 0.5 * zpsize;
    double xcenter = box.center().x();
    double ycenter = box.center().y();
    double zcenter = box.center().z();

    // determine index configuration: i->colums, j->rows
    int Ni = xd ? Nx : Ny;
    int Nj = zd ? Nz : Ny;

    // allocate result array with the appropriate size
    Array Tv(Ni * Nj);

    // calculate the results in parallel; perform only at the root processs
    auto parallel = probe->find<ParallelFactory>()->parallelRootOnly();
    parallel->call(Nj, [&Tv, ms, grid, units, xpsize, ypsize, zpsize, xbase, ybase, zbase, xd, yd, zd, xc, yc, zc,
                        Ni](size_t firstIndex, size_t numIndices) {
        for (size_t j = firstIndex; j != firstIndex + numIndices; ++j)
        {
            double z = zd ? (zbase + j * zpsize) : zc;
            for (int i = 0; i < Ni; i++)
            {
                double x = xd ? (xbase + i * xpsize) : xc;
                double y = yd ? (ybase + (zd ? i : j) * ypsize) : yc;
                int l = i + Ni * j;

                int m = grid->cellIndex(Position(x, y, z));
                if (m >= 0) Tv[l] = units->otemperature(ms->indicativeGasTemperature(m));
            }
        }
    });

    // get the name of the coordinate plane (xy, xz, or yz)
    string plane;
    if (xd) plane += "x";
    if (yd) plane += "y";
    if (zd) plane += "z";

    // write the results to a FITS file with an appropriate name
    string description = "gas temperatures in the " + plane + " plane";
    string filename = probe->itemName() + "_gas_T_" + plane;
    FITSInOut::write(probe, description, filename, Tv, units->utemperature(), Ni, Nj,
                     units->olength(xd ? xpsize : ypsize), units->olength(zd ? zpsize : ypsize),
                     units->olength(xd ? xcenter : ycenter), units->olength(zd ? zcenter : ycenter), units->ulength());
}

////////////////////////////////////////////////////////////////////

void PlanarGasTemperatureCutsProbe::probe()
{
    if (find<Configuration>()->hasMedium() && find<MediumSystem>()->hasGas())
    {
        writeGasTemperatureCut(this, 1, 1, 0, positionX(), positionY(), positionZ(), numPixelsX(), numPixelsY(),
                               numPixelsZ());
        writeGasTemperatureCut(this, 1, 0, 1, positionX(), positionY(), positionZ(), numPixelsX(), numPixelsY(),
                               numPixelsZ());
        writeGasTemperatureCut(this, 0, 1, 1, positionX(), positionY(), positionZ(), numPixelsX(), numPixelsY(),
                               numPixelsZ());
    }
}

////////////////////////////////////////////////////////////////////
