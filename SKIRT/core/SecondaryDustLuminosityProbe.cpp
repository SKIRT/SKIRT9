/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SecondaryDustLuminosityProbe.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "Indices.hpp"
#include "MediumSystem.hpp"
#include "ProbeFormBridge.hpp"
#include "SpecialFunctions.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

Probe::When SecondaryDustLuminosityProbe::when() const
{
    switch (probeAfter())
    {
        case ProbeAfter::Run: return When::Run;
        case ProbeAfter::Secondary: return When::Secondary;
    }
    return When::Run;
}

////////////////////////////////////////////////////////////////////

void SecondaryDustLuminosityProbe::probe()
{
    // locate the configuration and medium system; abort if there is no output for this probe
    auto config = find<Configuration>();
    if (!config->hasDustEmission()) return;
    auto ms = find<MediumSystem>();
    if (!ms->hasDust()) return;

    // get the units system and dust emission wavelength grid
    auto units = find<Units>();
    auto wlg = config->dustEmissionWLG();

    // construct lists of wavelength information in **output order**
    Array axis(wlg->numBins());  // representative wavelengths in output units
    Array cvol(wlg->numBins());  // conversion factors for luminosity volume density
    Array csrf(wlg->numBins());  // conversion factors for luminosity surface density (i.e. surface brightness)
    int outell = 0;
    for (int ell : Indices(wlg->numBins(), units->rwavelength()))
    {
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

    // get the extended dust emission wavelength grid including outer borders
    Array ewlg = wlg->extlambdav();

    // define the call-back function to retrieve a luminosity volume density spectrum in output ordering
    auto valueInCell = [ms, units, &ewlg, &cvol](int m) {
        // get the emmissivity spectrum for all dust in the cell with arbitrary normalization
        const Array& ev = ms->dustEmissionSpectrum(m);

        // calculate the normalization factor, assuming logarithmic interpolation (see NR::cdf2(...))
        double norm = 0.;
        size_t n = ewlg.size() - 1;
        for (size_t i = 0; i != n; ++i)
        {
            if (ev[i] > 0 && ev[i + 1] > 0)
            {
                double alpha = log(ev[i + 1] / ev[i]) / log(ewlg[i + 1] / ewlg[i]);
                norm += ev[i] * ewlg[i] * SpecialFunctions::gln(-alpha, ewlg[i + 1] / ewlg[i]);
            }
        }

        // determine the total dust luminosity volume density in the cell, adjusted with the normalization factor
        double front = norm > 0. ? ms->dustLuminosity(m) / ms->grid()->volume(m) / norm : 0.;

        // copy and renormalize the values in output order and units, omitting the outer borders
        int numBins = ewlg.size() - 2;
        Array Lv(numBins);
        int outell = 0;
        for (int ell : Indices(numBins, units->rwavelength()))
        {
            Lv[outell] = cvol[outell] * front * ev[ell + 1];  // skip left outer border
            outell++;
        }
        return Lv;
    };

    // construct a bridge and tell it to output the luminosity
    ProbeFormBridge bridge(this, form());
    bridge.writeQuantity("dust_L", "dust_S", units->umonluminosityvolumedensity(), units->usurfacebrightness(),
                         csrf[0] / cvol[0], "luminosity density", "surface brightness", axis, units->uwavelength(),
                         addColumnDefinitions, valueInCell);
}

////////////////////////////////////////////////////////////////////
