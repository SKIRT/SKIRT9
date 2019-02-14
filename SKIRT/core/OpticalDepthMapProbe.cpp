/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "OpticalDepthMapProbe.hpp"
#include "Array.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "FITSInOut.hpp"
#include "HomogeneousTransform.hpp"
#include "MaterialMix.hpp"
#include "Medium.hpp"
#include "MediumSystem.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "SpatialGridPath.hpp"
#include "StringUtils.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void OpticalDepthMapProbe::setupSelfBefore()
{
    Probe::setupSelfBefore();

    // verify consistency of vector properties
    Vec co(_Ox-_Cx, _Oy-_Cy, _Oz-_Cz); // vector in direction from crosshair to observer
    Vec up(_Ux, _Uy, _Uz);             // vector in the upward direction
    if (co.norm() < 1e-20) throw FATALERROR("Crosshair is too close to observer");
    if (up.norm() < 1e-20) throw FATALERROR("Upwards direction cannot be null vector");
    if (Vec::cross(co,up).norm() < 1e-20) throw FATALERROR("Upwards direction cannot be parallel to viewing direction");
}

////////////////////////////////////////////////////////////////////

namespace
{
    // private class to output the optical depth map for the given material type in parallel
    class WriteMap
    {
    private:
        // data members passed to constructor
        OpticalDepthMapProbe* probe;
        HomogeneousTransform& transform;
        MediumSystem* ms;
        MaterialMix::MaterialType type;
        string name;

        // data members derived in constructor
        int Nx, Ny;
        double lambda;

        // results -- resized in constructor
        Array tauv;

    public:
        // constructor
        WriteMap(OpticalDepthMapProbe* probe_, HomogeneousTransform& transform_, MediumSystem* ms_,
                 MaterialMix::MaterialType type_, string name_)
            : probe(probe_), transform(transform_), ms(ms_), type(type_), name(name_)
        {
            // set number of pixels using fixed aspect ratio
            Ny = probe->numPixelsY();
            Nx = 2*Ny;
            lambda = probe->wavelength();
            tauv.resize(Nx*Ny);
        }

        // the parallized loop body; calculates the results for a single line in the image
        void body(size_t firstIndex, size_t numIndices)
        {
            // construct a grid path with the user-configured starting point and an arbitary direction
            SpatialGridPath path(Position(probe->observerX(),probe->observerY(),probe->observerZ()), Direction());

            for (size_t j = firstIndex; j != firstIndex+numIndices; ++j)
            {
                double y = static_cast<double>(2*j+1)/static_cast<double>(Ny) - 1.;
                for (int i=0; i<Nx; i++)
                {
                    double x = static_cast<double>(2*i+1) / static_cast<double>(Nx) - 1.;

                    // convert viewport coordinates (-1 < x < 1  and -1 < y < 1 ) to spherical coordinates
                    double theta, phi;
                    bool inrange = probe->projection()->fromRectangleToSphere(x, y, theta, phi);

                    // if the deprojected direction is within range, compute the optical depth
                    if (inrange)
                    {
                        // set the path direction after transforming from observer to world coordinates
                        path.setDirection(Direction(transform.transform(Direction(theta, phi))));
                        tauv[i+Nx*j] = ms->opticalDepth(&path, lambda, type);
                    }
                }
            }
        }

        // write the results to a FITS file with an appropriate name
        void write()
        {
            Units* units = ms->find<Units>();
            string filename = probe->itemName() + "_" + name + "_tau";
            string description = "optical depth map at λ = "
                                  + StringUtils::toString(units->owavelength(probe->wavelength()))
                                  + " " + units->uwavelength();
            FITSInOut::write(probe, description, filename, tauv, "", Nx, Ny,
                             units->oposangle(2*M_PI/Nx), units->oposangle(M_PI/Ny), 0., 0.,
                             units->uposangle());
        }
    };

    // outputs the optical depth map for the given material type
    void writeMapForMaterialType(OpticalDepthMapProbe* probe, HomogeneousTransform& transform, MediumSystem* ms,
                                 MaterialMix::MaterialType type, string name)
    {
        // construct a private class instance to do the work (parallelized)
        WriteMap wm(probe, transform, ms, type, name);

        // perform the calculation at the root process only
        Parallel* parallel = probe->find<ParallelFactory>()->parallelRootOnly();
        parallel->call(probe->numPixelsY(), [&wm](size_t i ,size_t n) { wm.body(i, n); });

        // output the map
        wm.write();
    }
}

////////////////////////////////////////////////////////////////////

void OpticalDepthMapProbe::probeSetup()
{
    if (find<Configuration>()->hasMedium())
    {
        // locate the medium system
        auto ms = find<MediumSystem>();

        // setup the transformation from observer to world coordinates;
        // we need the inverse of the translation used by the all-sky instrument,
        // but without the translation because we will transform just the direction
        HomogeneousTransform transform;
        transform.rotateZ(0., 1.);
        transform.rotateX(0., -1.);

        Vec kn(_Ox-_Cx, _Oy-_Cy, _Oz-_Cz);
        kn /= kn.norm();
        double a = kn.x();
        double b = kn.y();
        double c = kn.z();

        double v = sqrt(b*b+c*c);
        if (v > 0.3)
        {
            double k = (b*b+c*c)*_Ux - a*b*_Uy - a*c*_Uz;
            double l = c*_Uy - b*_Uz;
            double u = sqrt(k*k+l*l);
            transform.rotateZ(l/u, k/u);
            transform.rotateY(v, a);
            transform.rotateX(c/v, b/v);
        }
        else
        {
            v = sqrt(a*a+c*c);
            double k = c*_Ux - a*_Uz;
            double l = (a*a+c*c)*_Uy - a*b*_Ux - b*c*_Uz;
            double u = sqrt(k*k+l*l);
            transform.rotateZ(l/u, k/u);
            transform.rotateX(v, b);
            transform.rotateY(c/v, a/v);
        }

        // write optical depth map for each material type, if present
        using Type = MaterialMix::MaterialType;
        if (ms->hasDust()) writeMapForMaterialType(this, transform, ms, Type::Dust, "dust");
        if (ms->hasElectrons()) writeMapForMaterialType(this, transform, ms, Type::Electrons, "elec");
        if (ms->hasGas()) writeMapForMaterialType(this, transform, ms, Type::Gas, "gas");
    }
}

////////////////////////////////////////////////////////////////////
