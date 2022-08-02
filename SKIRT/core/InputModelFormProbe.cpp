/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "InputModelFormProbe.hpp"
#include "Configuration.hpp"
#include "ImportedMedium.hpp"
#include "ImportedSource.hpp"
#include "MediumSystem.hpp"
#include "SourceSystem.hpp"

////////////////////////////////////////////////////////////////////

void InputModelFormProbe::probe()
{
    // handle sources together
    const auto& sources = find<SourceSystem>()->sources();
    if (!sources.empty())
    {
        // collect snapshots for all imported sources
        vector<const ImportedSource*> importedSources;
        vector<const Snapshot*> snapshots;
        for (auto source : sources)
        {
            auto importedSource = dynamic_cast<ImportedSource*>(source);
            if (!importedSource) break;
            importedSources.push_back(importedSource);
            snapshots.push_back(importedSource->snapshot());
        }

        // if all sources are imported, call subclass
        if (importedSources.size() == sources.size())
        {
            probeImportedSources(importedSources, snapshots);
        }
    }

    // loop over media
    if (find<Configuration>()->hasMedium())
    {
        const auto& media = find<MediumSystem>()->media();
        int numMedia = media.size();
        for (int h = 0; h != numMedia; ++h)
        {
            auto importedMedium = dynamic_cast<ImportedMedium*>(media[h]);
            if (importedMedium) probeImportedMedium(std::to_string(h), importedMedium, importedMedium->snapshot());
        }
    }
}

////////////////////////////////////////////////////////////////////

void InputModelFormProbe::probeImportedSources(const vector<const ImportedSource*>& /*sources*/,
                                               const vector<const Snapshot*>& /*snapshots*/)
{}

////////////////////////////////////////////////////////////////////

void InputModelFormProbe::probeImportedMedium(string /*sh*/, const ImportedMedium* /*medium*/,
                                              const Snapshot* /*snapshot*/)
{}

////////////////////////////////////////////////////////////////////
