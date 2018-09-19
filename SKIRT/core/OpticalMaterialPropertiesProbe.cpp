/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "OpticalMaterialPropertiesProbe.hpp"
#include "Configuration.hpp"
#include "InstrumentSystem.hpp"
#include "MaterialMix.hpp"
#include "Medium.hpp"
#include "MediumSystem.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // returns a human readable string describing the specified material type
    string materialForType(MaterialMix::MaterialType type)
    {
        switch(type)
        {
        case MaterialMix::MaterialType::Dust: return "dust";
        case MaterialMix::MaterialType::Electrons: return "electron";
        case MaterialMix::MaterialType::Gas: return "gas";
        }
        return string();  // to satisfy gcc compiler
    }

    // returns a human readable string describing a single entity in the specified material type
    string entityForType(MaterialMix::MaterialType type)
    {
        switch(type)
        {
        case MaterialMix::MaterialType::Dust: return "hydrogen atom";
        case MaterialMix::MaterialType::Electrons: return "electron";
        case MaterialMix::MaterialType::Gas: return "hydrogen atom";
        }
        return string();  // to satisfy gcc compiler
    }
}

////////////////////////////////////////////////////////////////////

void OpticalMaterialPropertiesProbe::probeSetup()
{
    if (find<Configuration>()->hasMedium())
    {
        auto units = find<Units>();

        // locate the medium system
        auto ms = find<MediumSystem>();
        int numMedia = ms->numMedia();

        // select "local" or default wavelength grid
        auto probeWaveGrid = wavelengthGrid() ? wavelengthGrid() : find<InstrumentSystem>()->find<WavelengthGrid>();
        int numWavelengths = probeWaveGrid->numBins();

        // create a seperate file for each medium
        for (int h=0; h!=numMedia; ++h)
        {
            // get the mix
            auto mix = ms->media()[h]->mix();

            // create a text file
            TextOutFile out(this, itemName() + "_opticalprops_" + std::to_string(h), "optical properties");

            // write the header
            out.writeLine("# Medium component " + std::to_string(h) + " -- " + materialForType(mix->materialType())
                          + " mass per " + entityForType(mix->materialType()) + ": "
                          + StringUtils::toString(units->obulkmass(mix->mass())) + " " + units->ubulkmass());
            out.addColumn("wavelength", units->uwavelength());
            out.addColumn("extinction cross section per " + entityForType(mix->materialType()), units->usection());
            out.addColumn("absorption cross section per " + entityForType(mix->materialType()), units->usection());
            out.addColumn("scattering cross section per " + entityForType(mix->materialType()), units->usection());
            out.addColumn("extinction mass coefficient", units->umasscoefficient());
            out.addColumn("absorption mass coefficient", units->umasscoefficient());
            out.addColumn("scattering mass coefficient", units->umasscoefficient());
            out.addColumn("scattering albedo");
            out.addColumn("scattering asymmetry parameter");

            // write the columns
            for (int ell=0; ell!=numWavelengths; ++ell)
            {
                double lambda = probeWaveGrid->wavelength(ell);

                out.writeRow(vector<double>({ units->owavelength(lambda),
                                              units->osection(mix->sectionExt(lambda)),
                                              units->osection(mix->sectionAbs(lambda)),
                                              units->osection(mix->sectionSca(lambda)),
                                              units->omasscoefficient(mix->sectionExt(lambda)/mix->mass()),
                                              units->omasscoefficient(mix->sectionAbs(lambda)/mix->mass()),
                                              units->omasscoefficient(mix->sectionSca(lambda)/mix->mass()),
                                              mix->albedo(lambda),
                                              mix->asymmpar(lambda) }));
            }
        }
    }
}

////////////////////////////////////////////////////////////////////
