/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "OpticalMaterialPropertiesProbe.hpp"
#include "Configuration.hpp"
#include "Indices.hpp"
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
        switch (type)
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
        switch (type)
        {
            case MaterialMix::MaterialType::Dust: return "hydrogen atom";
            case MaterialMix::MaterialType::Electrons: return "electron";
            case MaterialMix::MaterialType::Gas: return "hydrogen atom";
        }
        return string();  // to satisfy gcc compiler
    }
}

////////////////////////////////////////////////////////////////////

void OpticalMaterialPropertiesProbe::probe()
{
    if (find<Configuration>()->hasMedium())
    {
        auto units = find<Units>();

        // locate the medium system
        auto ms = find<MediumSystem>();
        int numMedia = ms->numMedia();

        // select "local" or default wavelength grid
        auto probeWavelengthGrid = find<Configuration>()->wavelengthGrid(wavelengthGrid());

        // create a seperate file for each medium
        for (int h = 0; h != numMedia; ++h)
        {
            // get the mix
            auto mix = ms->media()[h]->mix();

            // create a text file
            TextOutFile out(this, itemName() + "_opticalprops_" + std::to_string(h), "optical properties");

            // write the header
            out.writeLine("# Medium component " + std::to_string(h) + " -- " + materialForType(mix->materialType())
                          + " mass per " + entityForType(mix->materialType()) + ": "
                          + StringUtils::toString(units->obulkmass(mix->mass()), 'e', 9) + " " + units->ubulkmass());
            out.addColumn("wavelength; " + units->swavelength(), units->uwavelength());
            out.addColumn("extinction cross section per " + entityForType(mix->materialType()), units->usection());
            out.addColumn("absorption cross section per " + entityForType(mix->materialType()), units->usection());
            out.addColumn("scattering cross section per " + entityForType(mix->materialType()), units->usection());
            out.addColumn("extinction mass coefficient", units->umasscoefficient());
            out.addColumn("absorption mass coefficient", units->umasscoefficient());
            out.addColumn("scattering mass coefficient", units->umasscoefficient());
            out.addColumn("scattering albedo");
            out.addColumn("scattering asymmetry parameter");

            // write the columns
            for (int ell : Indices(probeWavelengthGrid->numBins(), units->rwavelength()))
            {
                double lambda = probeWavelengthGrid->wavelength(ell);
                double sigmaExt = mix->sectionExt(lambda);
                double sigmaAbs = mix->sectionAbs(lambda);
                double sigmaSca = mix->sectionSca(lambda);
                double albedo = sigmaExt ? sigmaSca / sigmaExt : 0.;
                out.writeRow(
                    vector<double>({units->owavelength(lambda), units->osection(sigmaExt), units->osection(sigmaAbs),
                                    units->osection(sigmaSca), units->omasscoefficient(sigmaExt / mix->mass()),
                                    units->omasscoefficient(sigmaAbs / mix->mass()),
                                    units->omasscoefficient(sigmaSca / mix->mass()), albedo, mix->asymmpar(lambda)}));
            }
        }
    }
}

////////////////////////////////////////////////////////////////////
