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
    // loop over sources
    const auto& sources = find<SourceSystem>()->sources();
    int numSources = sources.size();
    for (int h = 0; h != numSources; ++h)
    {
        auto source = dynamic_cast<ImportedSource*>(sources[h]);
        if (source) probeImportedSource(std::to_string(h), source, source->snapshot());
    }

    // loop over media
    if (find<Configuration>()->hasMedium())
    {
        const auto& media = find<MediumSystem>()->media();
        int numMedia = media.size();
        for (int h = 0; h != numMedia; ++h)
        {
            auto medium = dynamic_cast<ImportedMedium*>(media[h]);
            if (medium) probeImportedMedium(std::to_string(h), medium, medium->snapshot());
        }
    }
}

////////////////////////////////////////////////////////////////////

void InputModelFormProbe::probeImportedSource(string /*sh*/, const ImportedSource* /*source*/,
                                              const Snapshot* /*snapshot*/)
{}

////////////////////////////////////////////////////////////////////

void InputModelFormProbe::probeImportedMedium(string /*sh*/, const ImportedMedium* /*medium*/,
                                              const Snapshot* /*snapshot*/)
{}

////////////////////////////////////////////////////////////////////
