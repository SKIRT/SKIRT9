/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DefaultMediaDensityCutsProbe.hpp"
#include "Array.hpp"
#include "Box.hpp"
#include "FITSInOut.hpp"
#include "Log.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

// Private class to output FITS files with the theoretical mass density
// in one of the coordinate planes (xy, xz, or yz).
namespace
{
    // The image size in each direction, in pixels
    const int Np = 1024;

    class WriteDensity
    {
    private:
        // results -- sized to fit in constructor
        Array trhov;

        // data members initialized in constructor
        Probe* item;
        Geometry* geom;
        Units* units;
        Log* log;
        double xbase, ybase, zbase, xpsize, ypsize, zpsize, xcenter, ycenter, zcenter;

        // data members initialized in setup()
        bool xd, yd, zd; // direction of coordinate plane (110, 101, 011)
        string plane;    // name of the coordinate plane (xy, xz, yz)

    public:
        // constructor
        WriteDensity(Probe* item_, Geometry* geom_, const Box& bounds)
            : trhov(Np*Np), item(item_), geom(geom_), units(item_->find<Units>()), log(item_->find<Log>())
        {
             double xmin, ymin, zmin, xmax, ymax, zmax;
            bounds.extent(xmin,ymin,zmin,xmax,ymax,zmax);
            xpsize = (xmax-xmin)/Np;
            ypsize = (ymax-ymin)/Np;
            zpsize = (zmax-zmin)/Np;
            xbase = xmin + 0.5*xpsize;
            ybase = ymin + 0.5*ypsize;
            zbase = zmin + 0.5*zpsize;
            xcenter = (xmin+xmax)/2.0;
            ycenter = (ymin+ymax)/2.0;
            zcenter = (zmin+zmax)/2.0;
        }

        // setup for calculating a specific coordinate plane
        void setup(bool xdir, bool ydir, bool zdir)
        {
            xd = xdir;
            yd = ydir;
            zd = zdir;
            plane = "";
            if (xd) plane += "x";
            if (yd) plane += "y";
            if (zd) plane += "z";
            log->info(item->typeAndName() + " is calculating the media mass density in the " + plane + " plane...");
        }

        // the parallized loop body; calculates the results for a series of lines in the images
        void body(size_t firstIndex, size_t numIndices)
        {
            for (size_t j = firstIndex; j != firstIndex+numIndices; ++j)
            {
                double z = zd ? (zbase + j*zpsize) : 0.;
                for (int i=0; i<Np; i++)
                {
                    int l = i + Np*j;
                    double x = xd ? (xbase + i*xpsize) : 0.;
                    double y = yd ? (ybase + (zd ? i : j)*ypsize) : 0.;
                    Position bfr(x,y,z);
                    trhov[l] = units->omassvolumedensity(geom->density(bfr));
                }
            }
        }

        // write the results to two FITS files with appropriate names
        void write()
        {
            write(trhov, "theoretical density", item->itemName()+"_trho");
        }

    private:
        void write(const Array& rhov, string label, string prefix)
        {
            string filename = prefix + plane;
            FITSInOut::write(item, label + " in the " + plane + " plane", filename, rhov, Np, Np, 1,
                             units->olength(xd?xpsize:ypsize), units->olength(zd?zpsize:ypsize),
                             units->olength(xd?xcenter:ycenter), units->olength(zd?zcenter:ycenter),
                             units->umassvolumedensity(), units->ulength());
        }
    };
}

////////////////////////////////////////////////////////////////////

void DefaultMediaDensityCutsProbe::probeSetup()
{
    // construct a private class instance to do the work
    WriteDensity wd(this, _geometry, Box(_minX, _minY, _minZ, _maxX, _maxY, _maxZ));

    // configure parallelization; perform only at the root processs
    Parallel* parallel = find<ParallelFactory>()->parallelRootOnly();

    // get the dimension of the geometry
    int dimDust = _geometry->dimension();

    // For the xy plane (always)
    {
        wd.setup(1,1,0);
        parallel->call([&wd](size_t i ,size_t n) { wd.body(i, n); }, Np);
        wd.write();
    }

    // For the xz plane (only if dimension is at least 2)
    if (dimDust >= 2)
    {
        wd.setup(1,0,1);
        parallel->call([&wd](size_t i ,size_t n) { wd.body(i, n); }, Np);
        wd.write();
    }

    // For the yz plane (only if dimension is 3)
    if (dimDust == 3)
    {
        wd.setup(0,1,1);
        parallel->call([&wd](size_t i ,size_t n) { wd.body(i, n); }, Np);
        wd.write();
    }
}

////////////////////////////////////////////////////////////////////
