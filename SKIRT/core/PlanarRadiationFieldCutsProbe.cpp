/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PlanarRadiationFieldCutsProbe.hpp"
#include "Array.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "FITSInOut.hpp"
#include "InstrumentWavelengthGridProbe.hpp"
#include "MediumSystem.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void PlanarRadiationFieldCutsProbe::writeRadiationFieldCut(Probe* probe, bool xd, bool yd, bool zd, double xc,
                                                           double yc, double zc, int Nx, int Ny, int Nz)
{
    // locate relevant simulation items
    auto units = probe->find<Units>();
    auto ms = probe->find<MediumSystem>();
    auto grid = ms->grid();
    auto wavelengthGrid = probe->find<Configuration>()->radiationFieldWLG();

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

    // allocate result array with the appropriate size and initialize it to zero
    size_t size = Ni * Nj;
    Array Jvv(size * wavelengthGrid->numBins());

    // calculate the results in parallel; perform only at the root processs
    auto parallel = probe->find<ParallelFactory>()->parallelRootOnly();
    parallel->call(Nj, [&Jvv, units, ms, grid, wavelengthGrid, xpsize, ypsize, zpsize, xbase, ybase, zbase, xd, yd, zd,
                        xc, yc, zc, Ni, size](size_t firstIndex, size_t numIndices) {
        for (size_t j = firstIndex; j != firstIndex + numIndices; ++j)
        {
            double z = zd ? (zbase + j * zpsize) : zc;
            for (int i = 0; i < Ni; i++)
            {
                double x = xd ? (xbase + i * xpsize) : xc;
                double y = yd ? (ybase + (zd ? i : j) * ypsize) : yc;
                Position bfr(x, y, z);
                int m = grid->cellIndex(bfr);
                if (m >= 0)
                {
                    const Array& Jv = ms->meanIntensity(m);
                    for (int ell = 0; ell != wavelengthGrid->numBins(); ++ell)
                    {
                        size_t l = i + Ni * j + size * ell;
                        Jvv[l] = units->omeanintensity(wavelengthGrid->wavelength(ell), Jv[ell]);
                    }
                }
            }
        }
    });

    // get the name of the coordinate plane (xy, xz, or yz)
    string plane;
    if (xd) plane += "x";
    if (yd) plane += "y";
    if (zd) plane += "z";

    // write the results to a FITS file with an appropriate name
    string description = "mean intensity in the " + plane + " plane";
    string filename = probe->itemName() + "_J_" + plane;
    Array wavegrid(wavelengthGrid->numBins());
    for (int ell = 0; ell != wavelengthGrid->numBins(); ++ell)
        wavegrid[ell] = units->owavelength(wavelengthGrid->wavelength(ell));
    FITSInOut::write(probe, description, filename, Jvv, units->umeanintensity(), Ni, Nj,
                     units->olength(xd ? xpsize : ypsize), units->olength(zd ? zpsize : ypsize),
                     units->olength(xd ? xcenter : ycenter), units->olength(zd ? zcenter : ycenter), units->ulength(),
                     wavegrid, units->uwavelength());
}

////////////////////////////////////////////////////////////////////

void PlanarRadiationFieldCutsProbe::probeRun()
{
    if (find<Configuration>()->hasRadiationField())
    {
        writeRadiationFieldCut(this, 1, 1, 0, positionX(), positionY(), positionZ(), numPixelsX(), numPixelsY(),
                               numPixelsZ());
        writeRadiationFieldCut(this, 1, 0, 1, positionX(), positionY(), positionZ(), numPixelsX(), numPixelsY(),
                               numPixelsZ());
        writeRadiationFieldCut(this, 0, 1, 1, positionX(), positionY(), positionZ(), numPixelsX(), numPixelsY(),
                               numPixelsZ());

        // if requested, also output the wavelength grid
        if (writeWavelengthGrid())
        {
            InstrumentWavelengthGridProbe::writeWavelengthGrid(this, find<Configuration>()->radiationFieldWLG(),
                                                               itemName() + "_wavelengths",
                                                               "wavelengths for mean intensity");
        }
    }
}

////////////////////////////////////////////////////////////////////
