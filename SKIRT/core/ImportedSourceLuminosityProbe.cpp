/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedSourceLuminosityProbe.hpp"
#include "BandWavelengthGrid.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "EntityCollection.hpp"
#include "FatalError.hpp"
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

void ImportedSourceLuminosityProbe::probeImportedSources(const vector<const ImportedSource*>& sources,
                                                         const vector<const Snapshot*>& snapshots)
{
    // locate the units system
    auto units = find<Units>();

    // get the probe wavelength grid and determine the averaging style
    auto wlg = find<Configuration>()->wavelengthGrid(wavelengthGrid());
    auto dwlg = dynamic_cast<DisjointWavelengthGrid*>(wlg);
    auto bwlg = dynamic_cast<BandWavelengthGrid*>(wlg);
    enum class Style { Sample, Average, Convolve };
    Style style = Style::Sample;
    if (convolve() && !find<Configuration>()->oligochromatic())
    {
        if (dwlg)
            style = Style::Average;
        else if (bwlg)
            style = Style::Convolve;
        else
            throw FATALERROR("Convolve operation is not supported for wavelength grids of type " + wlg->type());
    }

    // construct lists of wavelength information in **output order**
    int numWaves = wlg->numBins();
    Array wave(numWaves);                // wavelengths in internal units (for style Sample)
    vector<Range> bin(numWaves);         // wavelength bin ranges in internal units (for style Average)
    vector<const Band*> band(numWaves);  // broadband for each wavelength bin (for style Convolve)
    Array axis(numWaves);                // wavelengths in output units
    Array cvol(numWaves);                // conversion factors for luminosity volume density
    Array csrf(numWaves);                // conversion factors for luminosity surface density (i.e. surface brightness)
    int outell = 0;
    for (int ell : Indices(numWaves, units->rwavelength()))
    {
        wave[outell] = wlg->wavelength(ell);
        if (style == Style::Average) bin[outell] = Range(wlg->leftBorder(ell), wlg->rightBorder(ell));
        if (style == Style::Convolve) band[outell] = bwlg->band(ell);
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

    // define a common call-back function to retrieve a compound luminosity value in one of two ways
    // depending on the value of the path argument (to avoid duplicating a lot of code):
    //  - a luminosity volume density value at a given position or along a given path
    //  - a surface brightness value along a given path
    auto valueAtPositionOrAlongPath = [&sources, &snapshots, numWaves, style, &wave, &bin, &band, &cvol,
                                       &csrf](bool path, Position bfr, Direction bfk) {
        // allocate an entity collection that can be reused for all queries in a given execution thread
        thread_local EntityCollection entities;

        // allocate array to store the result
        Array result(numWaves);

        // loop over all imported sources
        int numSources = sources.size();
        for (int h = 0; h != numSources; ++h)
        {
            // get the overlapping entities from the snapshot into the temporary collection
            if (path)
                snapshots[h]->getEntities(entities, bfr, bfk);
            else
                snapshots[h]->getEntities(entities, bfr);

            // loop over the entities
            for (const auto& entity : entities)
            {
                int m = entity.first;
                double w = entity.second;

                // loop over the wavelength bins
                for (int ell = 0; ell != numWaves; ++ell)
                {
                    // get the specific luminosity for each wavelength bin
                    double luminosity = 0.;
                    switch (style)
                    {
                        case Style::Sample: luminosity = sources[h]->specificLuminosity(wave[ell], m); break;
                        case Style::Average: luminosity = sources[h]->meanSpecificLuminosity(bin[ell], m); break;
                        case Style::Convolve: luminosity = sources[h]->meanSpecificLuminosity(band[ell], m); break;
                    }

                    // convert to the proper output value and accumulate into the result array
                    double c = path ? csrf[ell] : cvol[ell];
                    result[ell] += c * w * luminosity / snapshots[h]->volume(m);
                }
            }
        }
        return result;
    };

    // define the call-back function to retrieve a (compound) luminosity volume density value at a given position
    auto valueAtPosition = [valueAtPositionOrAlongPath](Position bfr) {
        return valueAtPositionOrAlongPath(false, bfr, Direction());
    };

    // define the call-back function to retrieve a surface brightness value along a given path
    auto valueAlongPath = [valueAtPositionOrAlongPath](Position bfr, Direction bfk) {
        return valueAtPositionOrAlongPath(true, bfr, bfk);
    };

    // construct a bridge and tell it to produce output
    ProbeFormBridge bridge(this, form());
    bridge.writeQuantity("L", "S", units->umonluminosityvolumedensity(), units->usurfacebrightness(),
                         "luminosity density", "surface brightness", axis, units->uwavelength(), addColumnDefinitions,
                         valueAtPosition, valueAlongPath);
}

////////////////////////////////////////////////////////////////////
