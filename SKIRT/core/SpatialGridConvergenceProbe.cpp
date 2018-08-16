/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpatialGridConvergenceProbe.hpp"
#include "Configuration.hpp"
#include "Medium.hpp"
#include "MediumSystem.hpp"
#include "SpatialGrid.hpp"
#include "SpatialGridPath.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // calculates the gridded optical depth at the given wavelength for the given material type
    // along the coordinate axis with the given direction
    double griddedOpticalDepth(MediumSystem* ms, double lambda, MaterialMix::MaterialType type, Direction axis)
    {
        double tau = 0;
        SpatialGridPath dgp(Position(0.,0.,0.), axis);
        ms->grid()->path(&dgp);
        tau += dgp.opticalDepth([ms,lambda,type](int m){ return ms->opacity(lambda, m, type); });
        dgp.setDirection(Direction(-axis.x(),-axis.y(),-axis.z()));
        ms->grid()->path(&dgp);
        tau += dgp.opticalDepth([ms,lambda,type](int m){ return ms->opacity(lambda, m, type); });
        return tau;
    }

    // writes two output lines with input and gridded values, respectively
    void writeValues(TextOutFile& out, double t, double g, double factor, string unit)
    {
        string tail;
        if (!unit.empty()) tail += " " + unit;
        out.writeLine("  Input:   " + StringUtils::toString(t*factor, 'g', 6) + tail);
        if (t > 0) tail += " (" + StringUtils::toString(100.*(g-t)/t,'f',2) + " %)";
        out.writeLine("  Gridded: " + StringUtils::toString(g*factor, 'g', 6) + tail);
    }

    // outputs the convergence info for the given material type
    void writeConvergenceForMaterialType(TextOutFile& out, MediumSystem* ms, double lambda,
                                         MaterialMix::MaterialType type, string name)
    {
        Units* units = ms->find<Units>();
        int numMedia = ms->numMedia();
        int numCells = ms->numCells();

        // output header
        out.writeLine("");
        out.writeLine("------ " + name + " ------");

        // calculate the true and gridded total mass
        double mass_t = 0.;
        double mass_g = 0.;
        for (int h=0; h!=numMedia; ++h)
        {
            if (ms->isMaterialType(type,h))
            {
                mass_t += ms->media()[h]->mass();
                for (int m=0; m!=numCells; ++m) mass_g += ms->massDensity(m,h)*ms->volume(m);
            }
        }

        // output the total mass
        out.writeLine("");
        out.writeLine("Total mass");
        writeValues(out, mass_t, mass_g, units->omass(1.), units->umass());

        // calculate the true optical depth along each of the coordinate axes
        double tau_X_t = 0.;
        double tau_Y_t = 0.;
        double tau_Z_t = 0.;
        for (int h=0; h!=numMedia; ++h)
        {
            if (ms->isMaterialType(type,h))
            {
                tau_X_t += ms->media()[h]->opticalDepthX(lambda);
                tau_Y_t += ms->media()[h]->opticalDepthY(lambda);
                tau_Z_t += ms->media()[h]->opticalDepthZ(lambda);
            }
        }

        // calculate the gridded optical depth along each of the coordinate axes
        double tau_X_g = griddedOpticalDepth(ms, lambda, type, Direction(1.,0.,0.));
        double tau_Y_g = griddedOpticalDepth(ms, lambda, type, Direction(0.,1.,0.));
        double tau_Z_g = griddedOpticalDepth(ms, lambda, type, Direction(0.,0.,1.));

        // output the optical depth along each of the coordinate axes
        out.writeLine("");
        out.writeLine("Optical depth at " + StringUtils::toString(units->owavelength(lambda), 'g')
                                + " " + units->uwavelength() + " along full X-axis");
        writeValues(out, tau_X_t, tau_X_g, 1., "");
        out.writeLine("");
        out.writeLine("Optical depth at " + StringUtils::toString(units->owavelength(lambda), 'g')
                                + " " + units->uwavelength() + " along full Y-axis");
        writeValues(out, tau_Y_t, tau_Y_g, 1., "");
        out.writeLine("");
        out.writeLine("Optical depth at " + StringUtils::toString(units->owavelength(lambda), 'g')
                                + " " + units->uwavelength() + " along full Z-axis");
        writeValues(out, tau_Z_t, tau_Z_g, 1., "");
    }
}

////////////////////////////////////////////////////////////////////

void SpatialGridConvergenceProbe::probeSetup()
{
    if (find<Configuration>()->hasMedia())
    {
        // locate the medium system
        auto ms = find<MediumSystem>();

        // create a text file
        TextOutFile out(this, itemName() + "_convergence", "spatial grid convergence info");
        out.writeLine("Spatial grid convergence information");
        out.writeLine("------------------------------------");

        // write info for each material type, if present
        using Type = MaterialMix::MaterialType;
        if (ms->hasDust()) writeConvergenceForMaterialType(out, ms, _wavelength, Type::Dust, "Dust");
        if (ms->hasElectrons()) writeConvergenceForMaterialType(out, ms, _wavelength, Type::Electrons, "Electrons");
        if (ms->hasGas()) writeConvergenceForMaterialType(out, ms, _wavelength, Type::Gas, "Gas");

        // write footer
        out.writeLine("");
        out.writeLine("------------------------------------");
    }
}

////////////////////////////////////////////////////////////////////
