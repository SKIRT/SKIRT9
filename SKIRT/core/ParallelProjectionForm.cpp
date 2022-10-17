/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ParallelProjectionForm.hpp"
#include "FITSInOut.hpp"
#include "Log.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "ProbeFormBridge.hpp"
#include "ProcessManager.hpp"
#include "Table.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void ParallelProjectionForm::writeQuantity(const ProbeFormBridge* bridge) const
{
    // get frame configuration
    int Nv = bridge->numValues();
    int Nxp = numPixelsX();
    int Nyp = numPixelsY();
    double xcent = centerX();
    double ycent = centerY();
    double xpmin = centerX() - 0.5 * fieldOfViewX();
    double xpsiz = fieldOfViewX() / numPixelsX();
    double ypmin = centerY() - 0.5 * fieldOfViewY();
    double ypsiz = fieldOfViewY() / numPixelsY();

    // calculate sines and cosines for the transformation from observer to model coordinates
    double costheta = cos(inclination());
    double sintheta = sin(inclination());
    double cosphi = cos(azimuth());
    double sinphi = sin(azimuth());
    double cosomega = cos(roll());
    double sinomega = sin(roll());

    // determine the observer axes as directions in model coordinates (inverse transform of frame instrument)
    // k_z is the direction from observer to model (i.e. left-handed coordinate frame)
    Direction kz(-cosphi * sintheta, -sinphi * sintheta, -costheta);
    Direction ky(-cosphi * costheta * cosomega - sinphi * sinomega, -sinphi * costheta * cosomega + cosphi * sinomega,
                 sintheta * cosomega);
    Direction kx(cosphi * costheta * sinomega - sinphi * cosomega, sinphi * costheta * sinomega + cosphi * cosomega,
                 -sintheta * sinomega);

    // get a distance that should be well outside of the model (but with a similar order of magnitude)
    double zp = 10. * (fieldOfViewX() + fieldOfViewY());

    // allocate result array with the appropriate size and initialize contents to zero
    Table<3> vvv(Nv, Nyp, Nxp);  // reverse index order to get proper data value ordering for FITSInOut::write()

    // calculate the results in parallel
    auto log = find<Log>();
    log->infoSetElapsed(Nyp);
    auto parallel = bridge->probe()->find<ParallelFactory>()->parallelDistributed();
    parallel->call(Nyp, [&vvv, bridge, xpmin, xpsiz, ypmin, ypsiz, kx, ky, kz, zp, costheta, sintheta, cosphi, sinphi,
                         cosomega, sinomega, Nxp, Nv, log](size_t firstIndex, size_t numIndices) {
        Array values(Nv);
        string progress = bridge->probe()->typeAndName() + " calculated projected pixels: ";

        // loop over pixels
        for (size_t j = firstIndex; j != firstIndex + numIndices; ++j)
        {
            for (int i = 0; i != Nxp; ++i)
            {
                // transform pixel indices to observer coordinates
                double xp = xpmin + (i + 0.5) * xpsiz;
                double yp = ypmin + (j + 0.5) * ypsiz;

                // transform observer coordinates to model coordinates (inverse transform of frame instrument)
                double xpp = sinomega * xp - cosomega * yp;
                double ypp = cosomega * xp + sinomega * yp;
                double zpp = zp;
                double x = cosphi * costheta * xpp - sinphi * ypp + cosphi * sintheta * zpp;
                double y = sinphi * costheta * xpp + cosphi * ypp + sinphi * sintheta * zpp;
                double z = -sintheta * xpp + costheta * zpp;

                // get the quantity value aggregated along a path leaving the pixel through the model
                bridge->valuesAlongPath(Position(x, y, z), kz, values);

                // store the result
                if (bridge->isVector())
                {
                    Vec v(values[0], values[1], values[2]);
                    vvv(0, j, i) = Vec::dot(v, kx);
                    vvv(1, j, i) = Vec::dot(v, ky);
                    vvv(2, j, i) = Vec::dot(v, kz);
                }
                else
                {
                    for (int p = 0; p != Nv; ++p) vvv(p, j, i) = values[p];
                }
            }
            log->infoIfElapsed(progress, 1);
        }
    });
    ProcessManager::sumToRoot(vvv.data(), true);

    // write the file
    auto units = bridge->units();
    FITSInOut::write(bridge->probe(), "parallel-projected " + bridge->projectedDescription(), bridge->projectedPrefix(),
                     vvv.data(), bridge->projectedUnit(), Nxp, Nyp, units->olength(xpsiz), units->olength(ypsiz),
                     units->olength(xcent), units->olength(ycent), units->ulength(), bridge->axis(),
                     bridge->axisUnit());
}

////////////////////////////////////////////////////////////////////
