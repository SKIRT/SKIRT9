/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedSourceLuminosityProbe.hpp"
#include "Configuration.hpp"
#include "EntityCollection.hpp"
#include "ImportedSource.hpp"
#include "Indices.hpp"
#include "ProbeFormBridge.hpp"
#include "Snapshot.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

Range ImportedSourceLuminosityProbe::wavelengthRange() const
{
    if (wavelengthGrid())
    {
        wavelengthGrid()->setup();
        return wavelengthGrid()->wavelengthRange();
    }
    return Range();
}

////////////////////////////////////////////////////////////////////

WavelengthGrid* ImportedSourceLuminosityProbe::materialWavelengthGrid() const
{
    return wavelengthGrid();
}

////////////////////////////////////////////////////////////////////

void ImportedSourceLuminosityProbe::probeImportedSource(string sh, const ImportedSource* source,
                                                        const Snapshot* snapshot)
{
    // locate the units system
    auto units = find<Units>();

    // get the probe wavelength grid
    auto wlg = find<Configuration>()->wavelengthGrid(wavelengthGrid());

    // construct the wavelength axis in output units and a list of corresponding unit conversion factors
    // ** all lists are in output order **
    int numWaves = wlg->numBins();
    Array wave(numWaves);  // wavelengths in internal units
    Array axis(numWaves);  // wavelengths in output units
    Array cvol(numWaves);  // conversion factors for luminosity volume density
    Array csrf(numWaves);  // conversion factors for luminosity surface density (i.e. surface brightness)
    int outell = 0;
    for (int ell : Indices(numWaves, units->rwavelength()))
    {
        wave[outell] = wlg->wavelength(ell);
        axis[outell] = units->owavelength(wlg->wavelength(ell));
        cvol[outell] = units->omonluminosityvolumedensity(wlg->wavelength(ell), 1.);
        csrf[outell] = units->osurfacebrightness(wlg->wavelength(ell), 1.) * (0.25 / M_PI);
        outell++;
    }

    // define the call-back function to add column definitions
    auto addColumnDefinitions = [axis, units](TextOutFile& outfile) {
        for (double outwave : axis)
        {
            // we assume that text files always contain values at a given position;
            // this will be incorrect if a new form would list projected quantities
            outfile.addColumn(units->smonluminosityvolumedensity() + " at " + units->swavelength() + " = "
                                  + StringUtils::toString(outwave, 'g') + " " + units->uwavelength(),
                              units->umonluminosityvolumedensity());
        }
    };

    // define the call-back function to retrieve a (compound) luminosity volume density value at a given position
    auto valueAtPosition = [source, snapshot, wave, numWaves, cvol](Position bfr) {
        thread_local EntityCollection entities;  // can be reused for all queries in a given execution thread
        snapshot->getEntities(entities, bfr);

        Array sumvw(numWaves);
        for (const auto& entity : entities)
        {
            int m = entity.first;
            double w = entity.second;
            for (int ell = 0; ell != numWaves; ++ell)
                sumvw[ell] += cvol[ell] * w * source->specificLuminosity(wave[ell], m) / snapshot->volume(m);
        }
        return sumvw;
    };

    // define the call-back function to retrieve a surface brightness value along a given path
    auto valueAlongPath = [source, snapshot, wave, numWaves, csrf](Position bfr, Direction bfk) {
        thread_local EntityCollection entities;  // can be reused for all queries in a given execution thread
        snapshot->getEntities(entities, bfr, bfk);

        Array sumvw(numWaves);
        for (const auto& entity : entities)
        {
            int m = entity.first;
            double w = entity.second;
            for (int ell = 0; ell != numWaves; ++ell)
                sumvw[ell] += csrf[ell] * w * source->specificLuminosity(wave[ell], m) / snapshot->volume(m);
        }
        return sumvw;
    };

    // construct a bridge and tell it to produce output
    ProbeFormBridge bridge(this, form());
    bridge.writeQuantity(sh + "_L", sh + "_S", units->umonluminosityvolumedensity(), units->usurfacebrightness(),
                         "luminosity density", "surface brightness", axis, units->uwavelength(), addColumnDefinitions,
                         valueAtPosition, valueAlongPath);
}

////////////////////////////////////////////////////////////////////
