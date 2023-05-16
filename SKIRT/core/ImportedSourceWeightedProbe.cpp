/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedSourceWeightedProbe.hpp"
#include "ImportedSource.hpp"
#include "Snapshot.hpp"
#include <unordered_map>

////////////////////////////////////////////////////////////////////

Range ImportedSourceWeightedProbe::wavelengthRange() const
{
    return weight() == Weight::Luminosity ? Range(wavelength(), wavelength()) : Range();
}

////////////////////////////////////////////////////////////////////

void ImportedSourceWeightedProbe::probeImportedSources(const vector<const ImportedSource*>& sources,
                                                       const vector<const Snapshot*>& snapshots)
{
    if (weight() == Weight::Luminosity)
    {
        // verify that all sources include the requested wavelength
        bool haveWavelength = true;
        for (auto source : sources) haveWavelength &= source->wavelengthRange().containsFuzzy(wavelength());
        if (haveWavelength)
        {
            // construct a map from snapshot to source
            std::unordered_map<const Snapshot*, const ImportedSource*> sourceForSnapshot;
            int numSources = sources.size();
            for (int h = 0; h != numSources; ++h) sourceForSnapshot.emplace(snapshots[h], sources[h]);

            // define function to retrieve weight
            auto getWeight = [sourceForSnapshot, this](const Snapshot* snapshot, int m) {
                return sourceForSnapshot.at(snapshot)->specificLuminosity(wavelength(), m) / snapshot->volume(m);
            };

            // perform probing
            probeImportedSourceWeighted("luminosity", snapshots, getWeight);
        }
    }

    if (weight() == Weight::InitialMass)
    {
        bool haveInitialMass = true;
        for (auto snapshot : snapshots) haveInitialMass &= snapshot->hasInitialMass();
        if (haveInitialMass)
        {
            auto getWeight = [](const Snapshot* snapshot, int m) {
                return snapshot->initialMass(m) / snapshot->volume(m);
            };
            probeImportedSourceWeighted("initial mass", snapshots, getWeight);
        }
    }

    if (weight() == Weight::CurrentMass)
    {
        bool haveCurrentMass = true;
        for (auto snapshot : snapshots) haveCurrentMass &= snapshot->hasCurrentMass();
        if (haveCurrentMass)
        {
            auto getWeight = [](const Snapshot* snapshot, int m) {
                return snapshot->currentMass(m) / snapshot->volume(m);
            };
            probeImportedSourceWeighted("current mass", snapshots, getWeight);
        }
    }
}

////////////////////////////////////////////////////////////////////
