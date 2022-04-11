/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "RadiationFieldProbe.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "Indices.hpp"
#include "InstrumentWavelengthGridProbe.hpp"
#include "MediumSystem.hpp"
#include "ProbeFormBridge.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

Probe::When RadiationFieldProbe::when() const
{
    return When::Run;
}

////////////////////////////////////////////////////////////////////

void RadiationFieldProbe::probe()
{
    if (find<Configuration>()->hasRadiationField())
    {
        // locate the medium system and units system
        auto ms = find<MediumSystem>();
        auto units = find<Units>();

        // get the radiation field wavelength grid
        auto wlg = find<Configuration>()->radiationFieldWLG();

        // construct the wavelength axis in output units and ordering
        Array axis(wlg->numBins());
        int outell = 0;
        for (int ell : Indices(wlg->numBins(), units->rwavelength()))
        {
            axis[outell++] = units->owavelength(wlg->wavelength(ell));
        }

        // define the call-back function to add column definitions
        auto addColumnDefinitions = [wlg, units](TextOutFile& outfile) {
            for (int ell : Indices(wlg->numBins(), units->rwavelength()))
            {
                outfile.addColumn(units->smeanintensity() + " at " + units->swavelength() + " = "
                                      + StringUtils::toString(units->owavelength(wlg->wavelength(ell)), 'g') + " "
                                      + units->uwavelength(),
                                  units->umeanintensity());
            }
        };

        // define the call-back function to retrieve a compound mean intensity value in output units and ordering
        auto valueInCell = [ms, wlg, units](int m) {
            const Array& Jv = ms->meanIntensity(m);
            Array outJv(wlg->numBins());
            int outell = 0;
            for (int ell : Indices(wlg->numBins(), units->rwavelength()))
            {
                outJv[outell++] = units->omeanintensity(wlg->wavelength(ell), Jv[ell]);
            }
            return outJv;
        };

        // define the call-back function to retrieve a cell weight; the units don't matter
        auto weightInCell = [ms](int m) { return ms->massDensity(m); };

        // construct a bridge and tell it to output the mean intensity
        ProbeFormBridge bridge(this, form());
        bridge.writeQuantity("J", units->umeanintensity(), "mean intensity", "density-weighted mean intensity", axis,
                             units->uwavelength(), addColumnDefinitions, valueInCell, weightInCell);

        // if requested, also output the wavelength grid
        if (writeWavelengthGrid())
        {
            InstrumentWavelengthGridProbe::writeWavelengthGrid(this, wlg, itemName() + "_wavelengths",
                                                               "wavelengths for mean intensity");
        }
    }
}

////////////////////////////////////////////////////////////////////
