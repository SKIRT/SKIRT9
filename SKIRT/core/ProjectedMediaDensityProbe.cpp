/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ProjectedMediaDensityProbe.hpp"
#include "Array.hpp"
#include "Configuration.hpp"
#include "FITSInOut.hpp"
#include "MediumSystem.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PathSegmentGenerator.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void ProjectedMediaDensityProbe::probeSetup()
{
    if (probeAfter() == ProbeAfter::Setup) probe();
}

////////////////////////////////////////////////////////////////////

void ProjectedMediaDensityProbe::probeRun()
{
    if (probeAfter() == ProbeAfter::Run) probe();
}

////////////////////////////////////////////////////////////////////

void ProjectedMediaDensityProbe::probe()
{
    if (find<Configuration>()->hasMedium())
    {
        // locate relevant simulation items
        auto units = find<Units>();
        auto ms = find<MediumSystem>();
        int numMedia = ms->numMedia();
        auto grid = ms->grid();

        // get frame configuration
        int Nxp = numPixelsX();
        int Nyp = numPixelsY();
        double xcent = centerX();
        double ycent = centerY();
        double xpmin = centerX() - 0.5 * fieldOfViewX();
        double xpsiz = fieldOfViewX() / numPixelsX();
        double ypmin = centerY() - 0.5 * fieldOfViewY();
        double ypsiz = fieldOfViewY() / numPixelsY();

        // determine direction from observer to model
        Direction bfkobs(inclination(), azimuth());
        bfkobs = -bfkobs;

        // calculate sines and cosines for the transformation from observer to model coordinates
        double costheta = cos(inclination());
        double sintheta = sin(inclination());
        double cosphi = cos(azimuth());
        double sinphi = sin(azimuth());
        double cosomega = cos(roll());
        double sinomega = sin(roll());

        // allocate result arrays with the appropriate size and initialize contents to zero
        Array dust_v;
        Array elec_v;
        Array gas_v;
        int size = Nxp * Nyp;
        if (ms->hasDust()) dust_v.resize(size);
        if (ms->hasElectrons()) elec_v.resize(size);
        if (ms->hasGas()) gas_v.resize(size);

        // calculate the results in parallel; perform only at the root processs
        auto parallel = find<ParallelFactory>()->parallelRootOnly();
        parallel->call(Nyp, [&dust_v, &elec_v, &gas_v, ms, numMedia, grid, xpmin, xpsiz, ypmin, ypsiz, bfkobs, costheta,
                             sintheta, cosphi, sinphi, cosomega, sinomega, Nxp](size_t firstIndex, size_t numIndices) {
            // allocate objects that can be reused for multiple pixels
            SpatialGridPath path;
            auto generator = grid->createPathSegmentGenerator();

            // loop over pixels
            for (size_t j = firstIndex; j != firstIndex + numIndices; ++j)
            {
                for (int i = 0; i != Nxp; ++i)
                {
                    int l = i + Nxp * j;

                    // transform pixel indices to observer coordinates
                    double xp = xpmin + (i + 0.5) * xpsiz;
                    double yp = ypmin + (j + 0.5) * ypsiz;
                    double zp = 2. * grid->boundingBox().diagonal();  // somewhere outside of the model

                    // transform observer coordinates to model coordinates (inverse transform of frame instrument)
                    double xpp = sinomega * xp - cosomega * yp;
                    double ypp = cosomega * xp + sinomega * yp;
                    double zpp = zp;
                    double x = cosphi * costheta * xpp - sinphi * ypp + cosphi * sintheta * zpp;
                    double y = sinphi * costheta * xpp + cosphi * ypp + sinphi * sintheta * zpp;
                    double z = -sintheta * xpp + costheta * zpp;

                    // initialize a path leaving the pixel through the model
                    path.setPosition(Position(x, y, z));
                    path.setDirection(bfkobs);

                    // calculate the column densities for each medium type along the path
                    generator->start(&path);
                    while (generator->next())
                    {
                        int m = generator->m();
                        if (m >= 0)
                        {
                            for (int h = 0; h != numMedia; ++h)
                            {
                                if (ms->isDust(h))
                                {
                                    dust_v[l] += ms->massDensity(m, h) * generator->ds();
                                }
                                else if (ms->isElectrons(h))
                                {
                                    elec_v[l] += ms->numberDensity(m, h) * generator->ds();
                                }
                                else if (ms->isGas(h))
                                {
                                    gas_v[l] += ms->numberDensity(m, h) * generator->ds();
                                }
                            }
                        }
                    }
                }
            }
        });

        // define a function to write a result array to a FITS file
        auto write = [this, units, xpsiz, ypsiz, xcent, ycent, Nxp, Nyp](Array& v, string label, string prefix,
                                                                         bool massDensity) {
            if (v.size())
            {
                // convert to output units
                v *= (massDensity ? units->omasssurfacedensity(1.) : units->onumbersurfacedensity(1.));
                string densityUnits = massDensity ? units->umasssurfacedensity() : units->unumbersurfacedensity();

                // write file
                string filename = itemName() + "_" + prefix;
                FITSInOut::write(this, label, filename, v, densityUnits, Nxp, Nyp, units->olength(xpsiz),
                                 units->olength(ypsiz), units->olength(xcent), units->olength(ycent), units->ulength());
            }
        };

        // write the results for each material type to a FITS file with an appropriate name
        write(dust_v, "dust surface mass density", "dust", true);
        write(elec_v, "electron column number density", "elec", false);
        write(gas_v, "gas column number density", "gas", false);
    }
}

////////////////////////////////////////////////////////////////////
