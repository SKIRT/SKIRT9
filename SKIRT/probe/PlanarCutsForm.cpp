/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PlanarCutsForm.hpp"
#include "Box.hpp"
#include "FITSInOut.hpp"
#include "FatalError.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "ProbeFormBridge.hpp"
#include "ProcessManager.hpp"
#include "Table.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void PlanarCutsForm::writeQuantity(const ProbeFormBridge* bridge) const
{
    Box box(minX(), minY(), minZ(), maxX(), maxY(), maxZ());
    writePlanarCut(bridge, 1, 1, 0, box, positionX(), positionY(), positionZ(), numPixelsX(), numPixelsY(),
                   numPixelsZ());
    writePlanarCut(bridge, 1, 0, 1, box, positionX(), positionY(), positionZ(), numPixelsX(), numPixelsY(),
                   numPixelsZ());
    writePlanarCut(bridge, 0, 1, 1, box, positionX(), positionY(), positionZ(), numPixelsX(), numPixelsY(),
                   numPixelsZ());
}

////////////////////////////////////////////////////////////////////

void PlanarCutsForm::writePlanarCut(const ProbeFormBridge* bridge, bool xd, bool yd, bool zd, const Box& box, double xp,
                                    double yp, double zp, int xn, int yn, int zn)
{
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

    // calculate the results in parallel
    auto parallel = bridge->probe()->find<ParallelFactory>()->parallelDistributed();
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
                if (bridge->isVector())
                {
                    Vec v(values[0], values[1], values[2]);
                    vvv(0, j, i) = xd ? v.x() : v.y();
                    vvv(1, j, i) = zd ? v.z() : v.y();
                    vvv(2, j, i) = xd ? (yd ? -v.z() : v.y()) : -v.x();
                }
                else
                {
                    for (int p = 0; p != Nv; ++p) vvv(p, j, i) = values[p];
                }
            }
        }
    });
    ProcessManager::sumToRoot(vvv.data(), true);

    // get the name of the coordinate plane (xy, xz, or yz)
    string plane;
    if (xd) plane += "x";
    if (yd) plane += "y";
    if (zd) plane += "z";

    // determine whether we are in the coordinate plane or parallel to it
    double offset = 0.;
    if (!xd) offset = xp;
    if (!yd) offset = yp;
    if (!zd) offset = yp;
    string position = offset ? " parallel to the " : " in the ";

    // write the file
    auto units = bridge->units();
    FITSInOut::write(bridge->probe(), bridge->description() + position + plane + " plane",
                     bridge->prefix() + "_" + plane, vvv.data(), bridge->unit(), Ni, Nj,
                     units->olength(xd ? xpsize : ypsize), units->olength(zd ? zpsize : ypsize),
                     units->olength(xd ? xcenter : ycenter), units->olength(zd ? zcenter : ycenter), units->ulength(),
                     bridge->axis(), bridge->axisUnit());
}

////////////////////////////////////////////////////////////////////
