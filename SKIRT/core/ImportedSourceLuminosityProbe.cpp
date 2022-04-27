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
    Array csrf(numWaves);  // conversion factors for luminosity surface density
    int outell = 0;
    for (int ell : Indices(numWaves, units->rwavelength()))
    {
        wave[outell] = wlg->wavelength(ell);
        axis[outell] = units->owavelength(wlg->wavelength(ell));
        cvol[outell] = 1.;  // TODO
        csrf[outell] = 1.;  // TODO
        outell++;
    }

    // define the call-back function to add column definitions
    auto addColumnDefinitions = [axis, units](TextOutFile& outfile) {
        for (double lambda : axis)
        {
            outfile.addColumn(units->smeanintensity() + " at " + units->swavelength() + " = "
                                  + StringUtils::toString(lambda, 'g') + " " + units->uwavelength(),
                              "W/m/m3");  // TODO
        }
    };

    // define a function to accumulate the luminosity spectrum for all entities in a collection
    auto accumulate = [source, snapshot, wave, numWaves](const EntityCollection& entities) {
        Array sumvw(numWaves);
        for (const auto& entity : entities)
        {
            int m = entity.first;
            double w = entity.second;
            for (int ell = 0; ell != numWaves; ++ell)
                sumvw[ell] += source->specificLuminosity(wave[ell], m) / snapshot->volume(m) * w;
        }
        return sumvw;
    };

    // define the call-back function to retrieve a luminosity volume density value at a given position
    auto valueAtPosition = [snapshot, accumulate](Position bfr) {
        thread_local EntityCollection entities;  // can be reused for all queries in a given execution thread
        snapshot->getEntities(entities, bfr);
        return accumulate(entities);
    };

    // define the call-back function to retrieve a luminosity surface density value along a given path
    auto valueAlongPath = [snapshot, accumulate](Position bfr, Direction bfk) {
        thread_local EntityCollection entities;  // can be reused for all queries in a given execution thread
        snapshot->getEntities(entities, bfr, bfk);
        return accumulate(entities);
    };

    // construct a bridge and tell it to produce output
    ProbeFormBridge bridge(this, form());
    bridge.writeQuantity(sh + "_L", sh + "_L", "W/m/m3", "W/m/m2", "luminosity volume density",
                         "luminosity surface density", axis, units->uwavelength(), addColumnDefinitions,
                         valueAtPosition, valueAlongPath);  // TODO
}

////////////////////////////////////////////////////////////////////
