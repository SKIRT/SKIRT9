/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AllSkyProjectionForm.hpp"
#include "FITSInOut.hpp"
#include "FatalError.hpp"
#include "HomogeneousTransform.hpp"
#include "Log.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "ProbeFormBridge.hpp"
#include "ProcessManager.hpp"
#include "Table.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void AllSkyProjectionForm::setupSelfBefore()
{
    GenericForm::setupSelfBefore();

    // verify consistency of vector properties
    Vec co(_Ox - _Cx, _Oy - _Cy, _Oz - _Cz);  // vector in direction from crosshair to observer
    Vec up(_Ux, _Uy, _Uz);                    // vector in the upward direction
    if (co.norm() < 1e-20) throw FATALERROR("Crosshair is too close to observer");
    if (up.norm() < 1e-20) throw FATALERROR("Upwards direction cannot be null vector");
    if (Vec::cross(co, up).norm() < 1e-20)
        throw FATALERROR("Upwards direction cannot be parallel to viewing direction");
}

////////////////////////////////////////////////////////////////////

void AllSkyProjectionForm::writeQuantity(const ProbeFormBridge* bridge) const
{
    // setup the transformation from observer to model coordinates;
    // we need the inverse of the transformation used by the all-sky instrument,
    // but without the translation because we will transform just the direction
    HomogeneousTransform transform;
    transform.rotateZ(0., -1.);
    transform.rotateX(0., 1.);

    Vec kn(_Ox - _Cx, _Oy - _Cy, _Oz - _Cz);
    kn /= kn.norm();
    double a = kn.x();
    double b = kn.y();
    double c = kn.z();

    double v = sqrt(b * b + c * c);
    if (v > 0.3)
    {
        double k = (b * b + c * c) * _Ux - a * b * _Uy - a * c * _Uz;
        double l = c * _Uy - b * _Uz;
        double u = sqrt(k * k + l * l);
        transform.rotateZ(l / u, k / u);
        transform.rotateY(v, a);
        transform.rotateX(c / v, b / v);
    }
    else
    {
        v = sqrt(a * a + c * c);
        double k = c * _Ux - a * _Uz;
        double l = (a * a + c * c) * _Uy - a * b * _Ux - b * c * _Uz;
        double u = sqrt(k * k + l * l);
        transform.rotateZ(l / u, k / u);
        transform.rotateX(v, b);
        transform.rotateY(c / v, a / v);
    }

    // get frame configuration using fixed aspect ratio
    int Nv = bridge->numValues();
    int Ny = numPixelsY();
    int Nx = 2 * Ny;

    // allocate result array with the appropriate size and initialize contents to zero
    Table<3> vvv(Nv, Ny, Nx);  // reverse index order to get proper data value ordering for FITSInOut::write()

    // calculate the results in parallel
    auto log = find<Log>();
    log->infoSetElapsed(Ny);
    auto parallel = bridge->probe()->find<ParallelFactory>()->parallelDistributed();
    parallel->call(Ny, [&vvv, this, bridge, &transform, Nx, Ny, Nv, log](size_t firstIndex, size_t numIndices) {
        Array values(Nv);
        string progress = bridge->probe()->typeAndName() + " calculated projected pixels: ";

        // loop over pixels
        for (size_t j = firstIndex; j != firstIndex + numIndices; ++j)
        {
            double y = static_cast<double>(2 * j + 1) / static_cast<double>(Ny) - 1.;
            for (int i = 0; i != Nx; ++i)
            {
                double x = static_cast<double>(2 * i + 1) / static_cast<double>(Nx) - 1.;

                // convert viewport coordinates (-1 < x < 1  and -1 < y < 1 ) to spherical coordinates
                double theta, phi;
                bool inrange = projection()->fromRectangleToSphere(x, y, theta, phi);

                // if the deprojected direction is within range, compute and store the projected quantity
                if (inrange)
                {
                    // transform the direction from observer to model coordinates
                    Direction kz = Direction(transform.transform(Direction(theta, phi)));

                    // get the projected quantity values
                    bridge->valuesAlongPath(Position(_Ox, _Oy, _Oz), kz, values);

                    // store the result
                    if (bridge->isVector())
                    {
                        // obtain the x and y axis directions in model coordinates
                        Vec ku(_Ux, _Uy, _Uz);
                        Vec kx = Vec::cross(ku, kz);
                        Vec ky = Vec::cross(kz, kx);
                        kx = kx / kx.norm();
                        ky = ky / ky.norm();

                        // project the vector on each of the axes
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
            }
            log->infoIfElapsed(progress, 1);
        }
    });
    ProcessManager::sumToRoot(vvv.data(), true);

    // write the file
    auto units = bridge->units();
    FITSInOut::write(bridge->probe(), "all-sky-projected " + bridge->projectedDescription(), bridge->projectedPrefix(),
                     vvv.data(), bridge->projectedUnit(), Nx, Ny, units->oposangle(-2 * M_PI / Nx),
                     units->oposangle(M_PI / Ny), 0., 0., units->uposangle(), bridge->axis(), bridge->axisUnit());
}

////////////////////////////////////////////////////////////////////
