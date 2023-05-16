/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SecondaryLineLuminosityProbe.hpp"
#include "Configuration.hpp"
#include "Indices.hpp"
#include "MediumSystem.hpp"
#include "ProbeFormBridge.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

Probe::When SecondaryLineLuminosityProbe::when() const
{
    switch (probeAfter())
    {
        case ProbeAfter::Run: return When::Run;
        case ProbeAfter::Secondary: return When::Secondary;
    }
    return When::Run;
}

////////////////////////////////////////////////////////////////////

void SecondaryLineLuminosityProbe::probe()
{
    if (find<Configuration>()->hasGasEmission())
    {
        // locate the medium system and units system
        auto ms = find<MediumSystem>();
        auto units = find<Units>();

        // loop over medium components with line emission
        for (int h : ms->gasMediumIndices())
        {
            auto mix = ms->mix(0, h);
            if (mix->hasLineEmission())
            {
                // get the central wavelengths
                const Array& centers = mix->lineEmissionCenters();

                // construct the wavelength axis in output units and ordering
                Array axis(centers.size());
                int outell = 0;
                for (int ell : Indices(centers, units->rwavelength()))
                {
                    axis[outell] = units->owavelength(centers[ell]);
                    outell++;
                }

                // define the call-back function to add column definitions
                auto addColumnDefinitions = [centers, units](TextOutFile& outfile) {
                    for (int ell : Indices(centers, units->rwavelength()))
                    {
                        outfile.addColumn("L/V for line at "
                                              + StringUtils::toString(units->owavelength(centers[ell]), 'g') + " "
                                              + units->uwavelength(),
                                          units->ubolluminosityvolumedensity());
                    }
                };

                // define the call-back function to retrieve a compound luminosity value in output ordering
                auto valueInCell = [ms, h, units](int m) {
                    Array Lv = ms->lineEmissionSpectrum(m, h) / ms->grid()->volume(m);
                    if (units->rwavelength()) std::reverse(begin(Lv), end(Lv));
                    return Lv;
                };

                // construct a bridge and tell it to output the luminosity
                ProbeFormBridge bridge(this, form());
                string sh = std::to_string(h);
                bridge.writeQuantity(sh + "_L", sh + "_S", "bolluminosityvolumedensity", "bolluminositysurfacedensity",
                                     "line luminosity volume density", "line luminosity surface density", axis,
                                     units->uwavelength(), addColumnDefinitions, valueInCell);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////
