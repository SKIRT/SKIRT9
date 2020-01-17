/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PlanarMediaDensityCutsProbe.hpp"
#include "Array.hpp"
#include "Configuration.hpp"
#include "FITSInOut.hpp"
#include "MediumSystem.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void PlanarMediaDensityCutsProbe::writeMediaDensityCuts(Probe* probe, bool xd, bool yd, bool zd, double xc, double yc,
                                                        double zc, int Nx, int Ny, int Nz)
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

    // allocate result arrays with the appropriate size and initialize contents to zero
    Array dust_tv, dust_gv;
    Array elec_tv, elec_gv;
    Array gas_tv, gas_gv;
    int size = Ni * Nj;
    if (ms->hasDust())
    {
        dust_tv.resize(size), dust_gv.resize(size);
    }
    if (ms->hasElectrons())
    {
        elec_tv.resize(size), elec_gv.resize(size);
    }
    if (ms->hasGas())
    {
        gas_tv.resize(size), gas_gv.resize(size);
    }

    // calculate the results in parallel; perform only at the root processs
    auto parallel = probe->find<ParallelFactory>()->parallelRootOnly();
    parallel->call(Nj, [&dust_tv, &dust_gv, &elec_tv, &elec_gv, &gas_tv, &gas_gv, ms, grid, xpsize, ypsize, zpsize,
                        xbase, ybase, zbase, xd, yd, zd, xc, yc, zc, Ni](size_t firstIndex, size_t numIndices) {
        int numMedia = ms->numMedia();
        for (size_t j = firstIndex; j != firstIndex + numIndices; ++j)
        {
            double z = zd ? (zbase + j * zpsize) : zc;
            for (int i = 0; i < Ni; i++)
            {
                int l = i + Ni * j;
                double x = xd ? (xbase + i * xpsize) : xc;
                double y = yd ? (ybase + (zd ? i : j) * ypsize) : yc;
                Position bfr(x, y, z);
                int m = grid->cellIndex(bfr);

                for (int h = 0; h != numMedia; ++h)
                {
                    if (ms->isDust(h))
                    {
                        dust_tv[l] += ms->media()[h]->massDensity(bfr);
                        if (m >= 0) dust_gv[l] += ms->massDensity(m, h);
                    }
                    else if (ms->isElectrons(h))
                    {
                        elec_tv[l] += ms->media()[h]->numberDensity(bfr);
                        if (m >= 0) elec_gv[l] += ms->numberDensity(m, h);
                    }
                    else if (ms->isGas(h))
                    {
                        gas_tv[l] += ms->media()[h]->numberDensity(bfr);
                        if (m >= 0) gas_gv[l] += ms->numberDensity(m, h);
                    }
                }
            }
        }
    });

    // define a function to write a result array to a FITS file
    auto write = [probe, units, xpsize, ypsize, zpsize, xcenter, ycenter, zcenter, xd, yd, zd, Ni,
                  Nj](Array& v, string label, string prefix, bool massDensity) {
        if (v.size())
        {
            // get the name of the coordinate plane (xy, xz, or yz)
            string plane;
            if (xd) plane += "x";
            if (yd) plane += "y";
            if (zd) plane += "z";

            // convert to output units
            v *= (massDensity ? units->omassvolumedensity(1.) : units->onumbervolumedensity(1.));
            string densityUnits = massDensity ? units->umassvolumedensity() : units->unumbervolumedensity();

            // write file
            string filename = probe->itemName() + "_" + prefix + "_" + plane;
            FITSInOut::write(probe, label + " in the " + plane + " plane", filename, v, densityUnits, Ni, Nj,
                             units->olength(xd ? xpsize : ypsize), units->olength(zd ? zpsize : ypsize),
                             units->olength(xd ? xcenter : ycenter), units->olength(zd ? zcenter : ycenter),
                             units->ulength());
        }
    };

    // write the results for each material type to two FITS files with appropriate names
    write(dust_tv, "dust theoretical density", "dust_t", true);
    write(dust_gv, "dust gridded density", "dust_g", true);
    write(elec_tv, "electron theoretical density", "elec_t", false);
    write(elec_gv, "electron gridded density", "elec_g", false);
    write(gas_tv, "gas theoretical density", "gas_t", false);
    write(gas_gv, "gas gridded density", "gas_g", false);
}

////////////////////////////////////////////////////////////////////

void PlanarMediaDensityCutsProbe::probeSetup()
{
    if (find<Configuration>()->hasMedium())
    {
        writeMediaDensityCuts(this, 1, 1, 0, positionX(), positionY(), positionZ(), numPixelsX(), numPixelsY(),
                              numPixelsZ());
        writeMediaDensityCuts(this, 1, 0, 1, positionX(), positionY(), positionZ(), numPixelsX(), numPixelsY(),
                              numPixelsZ());
        writeMediaDensityCuts(this, 0, 1, 1, positionX(), positionY(), positionZ(), numPixelsX(), numPixelsY(),
                              numPixelsZ());
    }
}

////////////////////////////////////////////////////////////////////
