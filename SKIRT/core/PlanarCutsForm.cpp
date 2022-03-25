/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PlanarCutsForm.hpp"
#include "FITSInOut.hpp"
#include "FatalError.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "ProbeFormBridge.hpp"
#include "SpatialGrid.hpp"
#include "Table.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void PlanarCutsForm::writeFile(const ProbeFormBridge* bridge) const
{
    writePlanarCut(bridge, 1, 1, 0, positionX(), positionY(), positionZ(), fieldOfViewX(), fieldOfViewY(),
                   fieldOfViewZ(), numPixelsX(), numPixelsY(), numPixelsZ());
    writePlanarCut(bridge, 1, 0, 1, positionX(), positionY(), positionZ(), fieldOfViewX(), fieldOfViewY(),
                   fieldOfViewZ(), numPixelsX(), numPixelsY(), numPixelsZ());
    writePlanarCut(bridge, 0, 1, 1, positionX(), positionY(), positionZ(), fieldOfViewX(), fieldOfViewY(),
                   fieldOfViewZ(), numPixelsX(), numPixelsY(), numPixelsZ());
}

////////////////////////////////////////////////////////////////////

void PlanarCutsForm::writePlanarCut(const ProbeFormBridge* bridge, bool xd, bool yd, bool zd, double xp, double yp,
                                    double zp, double xf, double yf, double zf, int xn, int yn, int zn)
{
    // get the spatial grid bounds if available, and verify that we have them if we need them
    Box box;
    if (bridge->grid())
        box = bridge->grid()->boundingBox();
    else if (!xf || !yf || !zf)
        throw FATALERROR("Field of view must be explicitly configured if the simulation has no media");

    // determine spatial domain
    double xmin = xf ? xp - xf / 2. : box.xmin();
    double ymin = yf ? yp - yf / 2. : box.ymin();
    double zmin = zf ? zp - zf / 2. : box.zmin();
    double xmax = xf ? xp + xf / 2. : box.xmax();
    double ymax = yf ? yp + yf / 2. : box.ymax();
    double zmax = zf ? zp + zf / 2. : box.zmax();
    box = Box(xmin, ymin, zmin, xmax, ymax, zmax);

    // determine spatial configuration
    double xpsize = box.xwidth() / xn;
    double ypsize = box.ywidth() / yn;
    double zpsize = box.zwidth() / zn;
    double xbase = box.xmin() + 0.5 * xpsize;
    double ybase = box.ymin() + 0.5 * ypsize;
    double zbase = box.zmin() + 0.5 * zpsize;
    double xcenter = box.center().x();
    double ycenter = box.center().y();
    double zcenter = box.center().z();

    // determine index configuration: i->colums, j->rows
    int Ni = xd ? xn : yn;
    int Nj = zd ? zn : yn;
    int Nv = bridge->numValues();

    // allocate result array with the appropriate size and initialize contents to zero
    Table<3> vvv(Nv, Nj, Ni);  // reverse index order to get proper data value ordering for FITSInOut::write()

    // calculate the results in parallel; perform only at the root process
    auto parallel = bridge->probe()->find<ParallelFactory>()->parallelRootOnly();
    parallel->call(Nj, [&vvv, bridge, xpsize, ypsize, zpsize, xbase, ybase, zbase, xd, yd, zd, xp, yp, zp, Ni,
                        Nv](size_t firstIndex, size_t numIndices) {
        Array values(Nv);
        for (size_t j = firstIndex; j != firstIndex + numIndices; ++j)
        {
            double z = zd ? (zbase + j * zpsize) : zp;
            for (int i = 0; i < Ni; i++)
            {
                double x = xd ? (xbase + i * xpsize) : xp;
                double y = yd ? (ybase + (zd ? i : j) * ypsize) : yp;
                bridge->valuesAtPosition(Position(x, y, z), values);
                for (int p = 0; p != Nv; ++p) vvv(p, j, i) = values[p];
            }
        }
    });

    // get the name of the coordinate plane (xy, xz, or yz)
    string plane;
    if (xd) plane += "x";
    if (yd) plane += "y";
    if (zd) plane += "z";

    // write file
    auto units = bridge->probe()->find<Units>();
    FITSInOut::write(bridge->probe(), bridge->description() + " in the " + plane + " plane", bridge->prefix() + "_" + plane,
                     vvv.data(),
                     unitstring, Ni, Nj, units->olength(xd ? xpsize : ypsize), units->olength(zd ? zpsize : ypsize),
                     units->olength(xd ? xcenter : ycenter), units->olength(zd ? zcenter : ycenter), units->ulength(),
                     NR::array(std::vector<double>({1., 2., 3.})), "1");
}

////////////////////////////////////////////////////////////////////
