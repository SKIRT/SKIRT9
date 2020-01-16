/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PlanarMagneticFieldCutsProbe.hpp"
#include "Configuration.hpp"
#include "FITSInOut.hpp"
#include "Medium.hpp"
#include "MediumSystem.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "Table.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void PlanarMagneticFieldCutsProbe::writeMagneticFieldCut(Probe* probe, bool xd, bool yd, bool zd, double xc, double yc,
                                                         double zc, int Nx, int Ny, int Nz)
{
    // get info on units
    auto units = probe->find<Units>();
    double unitfactor = units->omagneticfield(1.);
    string unitstring = units->umagneticfield();

    // locate the medium system
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

    // allocate result array with the appropriate size and initialize contents to zero
    Table<3> Bvv(3, Nj, Ni);  // reverse index order to get proper data value ordering for FITSInOut::write()

    // calculate the results in parallel; perform only at the root process
    auto parallel = probe->find<ParallelFactory>()->parallelRootOnly();
    parallel->call(Nj, [&Bvv, unitfactor, ms, grid, xpsize, ypsize, zpsize, xbase, ybase, zbase, xd, yd, zd, xc, yc, zc,
                        Ni](size_t firstIndex, size_t numIndices) {
        for (size_t j = firstIndex; j != firstIndex + numIndices; ++j)
        {
            double z = zd ? (zbase + j * zpsize) : zc;
            for (int i = 0; i < Ni; i++)
            {
                double x = xd ? (xbase + i * xpsize) : xc;
                double y = yd ? (ybase + (zd ? i : j) * ypsize) : yc;

                int m = grid->cellIndex(Position(x, y, z));
                if (m >= 0)
                {
                    Vec B = ms->magneticField(m);
                    Bvv(0, j, i) = (xd ? B.x() : B.y()) * unitfactor;
                    Bvv(1, j, i) = (zd ? B.z() : B.y()) * unitfactor;
                    Bvv(2, j, i) = (xd ? (yd ? -B.z() : B.y()) : -B.x()) * unitfactor;
                }
            }
        }
    });

    // get the name of the coordinate plane (xy, xz, or yz)
    string plane;
    if (xd) plane += "x";
    if (yd) plane += "y";
    if (zd) plane += "z";

    // write file
    FITSInOut::write(probe, "magnetic field in the " + plane + " plane", probe->itemName() + "_B_" + plane, Bvv.data(),
                     unitstring, Ni, Nj, units->olength(xd ? xpsize : ypsize), units->olength(zd ? zpsize : ypsize),
                     units->olength(xd ? xcenter : ycenter), units->olength(zd ? zcenter : ycenter), units->ulength(),
                     NR::array(std::vector<double>({1., 2., 3.})), "1");
}

////////////////////////////////////////////////////////////////////

void PlanarMagneticFieldCutsProbe::probeSetup()
{
    if (find<Configuration>()->hasMagneticField())
    {
        writeMagneticFieldCut(this, 1, 1, 0, positionX(), positionY(), positionZ(), numPixelsX(), numPixelsY(),
                              numPixelsZ());
        writeMagneticFieldCut(this, 1, 0, 1, positionX(), positionY(), positionZ(), numPixelsX(), numPixelsY(),
                              numPixelsZ());
        writeMagneticFieldCut(this, 0, 1, 1, positionX(), positionY(), positionZ(), numPixelsX(), numPixelsY(),
                              numPixelsZ());
    }
}

////////////////////////////////////////////////////////////////////
