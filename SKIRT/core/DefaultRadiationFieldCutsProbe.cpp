/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DefaultRadiationFieldCutsProbe.hpp"
#include "Array.hpp"
#include "Configuration.hpp"
#include "FITSInOut.hpp"
#include "Medium.hpp"
#include "MediumSystem.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "SpatialGrid.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

// Private class to output a FITS file with the mean intensity of the radiation field
// in each of the coordinate planes (xy, xz, or yz).
namespace
{
    // The image size in each direction, in pixels
    const int Np = 1024;

    class WriteMeanIntensityCut
    {
    private:
        // data members initialized in constructor
        Probe* item;
        MediumSystem* ms;
        SpatialGrid* grid;
        WavelengthGrid* wavelengthGrid;
        Units* units;
        double xbase, ybase, zbase, xpsize, ypsize, zpsize, xcenter, ycenter, zcenter;

        // data members initialized in setup()
        bool xd, yd, zd;  // direction of coordinate plane (110, 101, 011)
        string plane;     // name of the coordinate plane (xy, xz, yz)

        // results vector, properly sized in constructor and initialized to zero in setup()
        Array Jvv;

    public:
        // constructor
        WriteMeanIntensityCut(Probe* item_, MediumSystem* ms_)
            : item(item_), ms(ms_), grid(ms_->grid()),
              wavelengthGrid(ms_->find<Configuration>()->radiationFieldWavelengthGrid()), units(ms_->find<Units>())
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

            Jvv.resize(Np*Np*wavelengthGrid->numBins());
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

            Jvv = 0.;  // initialize all values to zero to facilitate the code in body()
        }

        // the parallized loop body; calculates the results for a series of lines in the images
        void body(size_t firstIndex, size_t numIndices)
        {
            for (size_t j = firstIndex; j != firstIndex+numIndices; ++j)
            {
                double z = zd ? (zbase + j*zpsize) : 0.;
                for (int i=0; i<Np; i++)
                {
                    double x = xd ? (xbase + i*xpsize) : 0.;
                    double y = yd ? (ybase + (zd ? i : j)*ypsize) : 0.;
                    Position bfr(x,y,z);
                    int m = grid->cellIndex(bfr);
                    if (m>=0)
                    {
                        const Array& Jv = ms->meanIntensity(m);
                        for (int ell=0; ell!=wavelengthGrid->numBins(); ++ell)
                        {
                            int l = i + Np*j + Np*Np*ell;
                            Jvv[l] = units->osurfacebrightnessWavelength(wavelengthGrid->wavelength(ell), Jv[ell]);
                        }
                    }
                }
            }
        }

        // Write the results to a FITS file with an appropriate name
        void write()
        {
            string description = "mean intensity in the " + plane + " plane";
            string filename = item->itemName() + "_J_" + plane;
            FITSInOut::write(item, description, filename, Jvv, Np, Np, wavelengthGrid->numBins(),
                             units->olength(xd?xpsize:ypsize), units->olength(zd?zpsize:ypsize),
                             units->olength(xd?xcenter:ycenter), units->olength(zd?zcenter:ycenter),
                             units->usurfacebrightness(), units->ulength());
        }
    };
}

////////////////////////////////////////////////////////////////////

void DefaultRadiationFieldCutsProbe::probeRun()
{
    if (find<Configuration>()->hasRadiationField())
    {
        // locate the medium system
        auto ms = find<MediumSystem>();

        // configure parallelization; perform only at the root processs
        Parallel* parallel = find<ParallelFactory>()->parallelRootOnly();

        // construct a private class instance to do the work
        WriteMeanIntensityCut wmi(this, ms);

        // for the xy plane (always)
        {
            wmi.setup(1,1,0);
            parallel->call(Np, [&wmi](size_t i ,size_t n) { wmi.body(i, n); });
            wmi.write();
        }

        // for the xz plane (only if dimension is at least 2)
        int dimension = ms->dimension();
        if (dimension >= 2)
        {
            wmi.setup(1,0,1);
            parallel->call(Np, [&wmi](size_t i ,size_t n) { wmi.body(i, n); });
            wmi.write();
        }

        // for the yz plane (only if dimension is 3)
        if (dimension == 3)
        {
            wmi.setup(0,1,1);
            parallel->call(Np, [&wmi](size_t i ,size_t n) { wmi.body(i, n); });
            wmi.write();
        }
    }
}

////////////////////////////////////////////////////////////////////
