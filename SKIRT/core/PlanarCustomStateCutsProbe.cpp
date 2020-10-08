/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PlanarCustomStateCutsProbe.hpp"
#include "Array.hpp"
#include "Configuration.hpp"
#include "FITSInOut.hpp"
#include "MediumSystem.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void PlanarCustomStateCutsProbe::writeCustomStateCuts(Probe* probe, bool xd, bool yd, bool zd, double xc, double yc,
                                                      double zc, int Nx, int Ny, int Nz)
{
    // locate relevant simulation items
    auto units = probe->find<Units>();
    auto ms = probe->find<MediumSystem>();
    auto grid = ms->grid();
    int numMedia = ms->numMedia();

    // get the custom state variable descriptors for each medium and abort if none of the media have custom state
    vector<vector<StateVariable>> descriptors(numMedia);
    bool hasCustomState = false;
    for (int h = 0; h != numMedia; ++h)
    {
        auto mix = ms->media()[h]->mix();
        for (auto candidate : mix->specificStateVariableInfo())
        {
            if (candidate.identifier() == StateVariable::Identifier::Custom)
            {
                descriptors[h].push_back(candidate);
                hasCustomState = true;
            }
        }
    }
    if (!hasCustomState) return;

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
    vector<Array> results(numMedia);
    size_t size = Ni * Nj;
    for (int h = 0; h != numMedia; ++h) results[h].resize(descriptors[h].size() * size);

    // calculate the results in parallel; perform only at the root processs
    auto parallel = probe->find<ParallelFactory>()->parallelRootOnly();
    parallel->call(Nj, [&results, &descriptors, ms, numMedia, grid, xpsize, ypsize, zpsize, xbase, ybase, zbase, xd, yd,
                        zd, xc, yc, zc, Ni, size](size_t firstIndex, size_t numIndices) {
        for (size_t j = firstIndex; j != firstIndex + numIndices; ++j)
        {
            double z = zd ? (zbase + j * zpsize) : zc;
            for (int i = 0; i < Ni; i++)
            {
                double x = xd ? (xbase + i * xpsize) : xc;
                double y = yd ? (ybase + (zd ? i : j) * ypsize) : yc;
                Position bfr(x, y, z);
                int m = grid->cellIndex(bfr);

                for (int h = 0; h != numMedia; ++h)
                {
                    int numCustom = descriptors[h].size();
                    for (int c = 0; c != numCustom; ++c)
                    {
                        size_t l = i + Ni * j + size * c;
                        results[h][l] = ms->custom(m, h, descriptors[h][c].customIndex());
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

    // write the results to FITS files with appropriate names
    for (int h = 0; h != numMedia; ++h)
    {
        int numCustom = descriptors[h].size();
        if (numCustom)
        {
            string description = "custom state variables in the " + plane + " plane for medium " + std::to_string(h);
            string filename = probe->itemName() + "_customstate_" + plane + "_" + std::to_string(h);
            Array cgrid(numCustom);
            for (int c = 0; c != numCustom; ++c) cgrid[c] = descriptors[h][c].customIndex();
            FITSInOut::write(probe, description, filename, results[h], units->umeanintensity(), Ni, Nj,
                             units->olength(xd ? xpsize : ypsize), units->olength(zd ? zpsize : ypsize),
                             units->olength(xd ? xcenter : ycenter), units->olength(zd ? zcenter : ycenter),
                             units->ulength(), cgrid, "indices");
        }
    }
}

////////////////////////////////////////////////////////////////////

void PlanarCustomStateCutsProbe::probeSetup()
{
    if (probeAfter() == ProbeAfter::Setup) probe();
}

////////////////////////////////////////////////////////////////////

void PlanarCustomStateCutsProbe::probeRun()
{
    if (probeAfter() == ProbeAfter::Run) probe();
}

////////////////////////////////////////////////////////////////////

void PlanarCustomStateCutsProbe::probe()
{
    if (find<Configuration>()->hasMedium())
    {
        writeCustomStateCuts(this, 1, 1, 0, positionX(), positionY(), positionZ(), numPixelsX(), numPixelsY(),
                             numPixelsZ());
        writeCustomStateCuts(this, 1, 0, 1, positionX(), positionY(), positionZ(), numPixelsX(), numPixelsY(),
                             numPixelsZ());
        writeCustomStateCuts(this, 0, 1, 1, positionX(), positionY(), positionZ(), numPixelsX(), numPixelsY(),
                             numPixelsZ());
    }
}

////////////////////////////////////////////////////////////////////
