/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedSourceWeightedProbe.hpp"
#include "ImportedSource.hpp"
#include "Snapshot.hpp"

////////////////////////////////////////////////////////////////////

Range ImportedSourceWeightedProbe::wavelengthRange() const
{
    return weight() == Weight::Luminosity ? Range(wavelength(), wavelength()) : Range();
}

////////////////////////////////////////////////////////////////////

void ImportedSourceWeightedProbe::probeImportedSource(string sh, const ImportedSource* source, const Snapshot* snapshot)
{
    if (weight() == Weight::Luminosity && source->wavelengthRange().containsFuzzy(wavelength()))
    {
        probeImportedSourceWeighted(sh, "luminosity", snapshot, [source, snapshot, this](int m) {
            return source->specificLuminosity(wavelength(), m) / snapshot->volume(m);
        });
    }
    if (weight() == Weight::InitialMass && snapshot->hasInitialMass())
    {
        probeImportedSourceWeighted(sh, "initial mass", snapshot,
                                    [snapshot](int m) { return snapshot->initialMass(m) / snapshot->volume(m); });
    }
    if (weight() == Weight::CurrentMass && snapshot->hasCurrentMass())
    {
        probeImportedSourceWeighted(sh, "current mass", snapshot,
                                    [snapshot](int m) { return snapshot->currentMass(m) / snapshot->volume(m); });
    }
}

////////////////////////////////////////////////////////////////////
