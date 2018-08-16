/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DefaultMediaDensityCutsProbe.hpp"
#include "Array.hpp"
#include "Configuration.hpp"
#include "FITSInOut.hpp"
#include "Log.hpp"
#include "Medium.hpp"
#include "MediumSystem.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "SpatialGrid.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

// Private class to output FITS files with the theoretical and grid density
// for each material type in one of the coordinate planes (xy, xz, or yz).
namespace
{
    // The image size in each direction, in pixels
    const int Np = 1024;

    class WriteDensity
    {
    private:
        // data members initialized in constructor
        Probe* item;
        MediumSystem* ms;
        SpatialGrid* grid;
        Units* units;
        Log* log;
        double xbase, ybase, zbase, xpsize, ypsize, zpsize, xcenter, ycenter, zcenter;

        // data members initialized in setup()
        bool xd, yd, zd; // direction of coordinate plane (110, 101, 011)
        string plane;    // name of the coordinate plane (xy, xz, yz)

        // results -- resized and initialized to zero in setup()
        Array dust_tv, dust_gv;
        Array elec_tv, elec_gv;
        Array gas_tv, gas_gv;

    public:
        // constructor
        WriteDensity(Probe* item_, MediumSystem* ms_)
            : item(item_), ms(ms_), grid(ms_->grid()),
              units(ms_->find<Units>()), log(ms_->find<Log>())
        {
            double xmin, ymin, zmin, xmax, ymax, zmax;
            grid->boundingBox().extent(xmin,ymin,zmin,xmax,ymax,zmax);
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
            log->info(item->typeAndName() + " is calculating the media density in the " + plane + " plane...");

            if (ms->hasDust()) { dust_tv.resize(Np*Np), dust_gv.resize(Np*Np); }
            if (ms->hasElectrons()) { elec_tv.resize(Np*Np), elec_gv.resize(Np*Np); }
            if (ms->hasGas()) { gas_tv.resize(Np*Np), gas_gv.resize(Np*Np); }
        }

        // the parallized loop body; calculates the results for a series of lines in the images
        void body(size_t firstIndex, size_t numIndices)
        {
            int numMedia = ms->numMedia();
            for (size_t j = firstIndex; j != firstIndex+numIndices; ++j)
            {
                double z = zd ? (zbase + j*zpsize) : 0.;
                for (int i=0; i<Np; i++)
                {
                    int l = i + Np*j;
                    double x = xd ? (xbase + i*xpsize) : 0.;
                    double y = yd ? (ybase + (zd ? i : j)*ypsize) : 0.;
                    Position bfr(x,y,z);
                    int m = grid->cellIndex(bfr);

                    for (int h=0; h!=numMedia; ++h)
                    {
                        if (ms->isDust(h))
                        {
                            dust_tv[l] += ms->media()[h]->massDensity(bfr);
                            if (m>=0) dust_gv[l] += ms->massDensity(m,h);
                        }
                        else if (ms->isElectrons(h))
                        {
                            elec_tv[l] += ms->media()[h]->numberDensity(bfr);
                            if (m>=0) elec_gv[l] += ms->numberDensity(m,h);
                        }
                        else if (ms->isGas(h))
                        {
                            gas_tv[l] += ms->media()[h]->numberDensity(bfr);
                            if (m>=0) gas_gv[l] += ms->numberDensity(m,h);
                        }
                    }
                }
            }
        }

        // write the results to two FITS files with appropriate names
        void write()
        {
            write(dust_tv, "dust theoretical density", "dust_t", true);
            write(dust_gv, "dust gridded density", "dust_g", true);
            write(elec_tv, "electron theoretical density", "elec_t", false);
            write(elec_gv, "electron gridded density", "elec_g", false);
            write(gas_tv, "gas theoretical density", "gas_t", false);
            write(gas_gv, "gas gridded density", "gas_g", false);
        }

    private:
        void write(Array& v, string label, string prefix, bool massDensity)
        {
            if (v.size())
            {
                // convert to output units
                v *= (massDensity ? units->omassvolumedensity(1.) : units->onumbervolumedensity(1.));
                string densityUnits = massDensity ? units->umassvolumedensity() : units->unumbervolumedensity();

                // write file
                string filename = item->itemName() + "_" + prefix + "_" + plane;
                FITSInOut::write(item, label + " in the " + plane + " plane", filename, v, Np, Np, 1,
                                 units->olength(xd?xpsize:ypsize), units->olength(zd?zpsize:ypsize),
                                 units->olength(xd?xcenter:ycenter), units->olength(zd?zcenter:ycenter),
                                 densityUnits, units->ulength());
            }
        }
    };
}

////////////////////////////////////////////////////////////////////

void DefaultMediaDensityCutsProbe::probeSetup()
{
    if (find<Configuration>()->hasMedia())
    {
        // locate the medium system
        auto ms = find<MediumSystem>();

        // configure parallelization; perform only at the root processs
        Parallel* parallel = find<ParallelFactory>()->parallelRootOnly();

        // construct a private class instance to do the work
        WriteDensity wd(this, ms);

        // for the xy plane (always)
        {
            wd.setup(1,1,0);
            parallel->call(Np, [&wd](size_t i ,size_t n) { wd.body(i, n); });
            wd.write();
        }

        // for the xz plane (only if dimension is at least 2)
        int dimension = ms->dimension();
        if (dimension >= 2)
        {
            wd.setup(1,0,1);
            parallel->call(Np, [&wd](size_t i ,size_t n) { wd.body(i, n); });
            wd.write();
        }

        // for the yz plane (only if dimension is 3)
        if (dimension == 3)
        {
            wd.setup(0,1,1);
            parallel->call(Np, [&wd](size_t i ,size_t n) { wd.body(i, n); });
            wd.write();
        }
    }
}

////////////////////////////////////////////////////////////////////
