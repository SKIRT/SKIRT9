/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "InputModelProbe.hpp"
#include "Configuration.hpp"
#include "EntityCollection.hpp"
#include "ImportedMedium.hpp"
#include "ImportedSource.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "Snapshot.hpp"
#include "SourceSystem.hpp"
#include "StringUtils.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    void probeSource(const SimulationItem* item, int index, const ImportedSource* source)
    {
        item->find<Log>()->info("== " + source->type() + " " + std::to_string(index) + " ==");
        auto snapshot = source->snapshot();

        EntityCollection entities;
        snapshot->getEntities(entities, Position(), Direction(1,0,0));

        if (snapshot->hasTemperature())
        {
            double sum = 0.;
            double sumw = 0.;
            for (const auto& entity : entities)
            {
                double value = snapshot->temperature(entity.first);
                double weight = entity.second;
                sum += value * weight;
                sumw += weight;
            }
            if (sumw) sum /= sumw;

            item->find<Log>()->info("   Temperature: " + StringUtils::toString(sum) + " K");
        }
        if (snapshot->hasMetallicity())
        {
            double sum = 0.;
            double sumw = 0.;
            for (const auto& entity : entities)
            {
                double value = snapshot->metallicity(entity.first);
                double weight = entity.second;
                sum += value * weight;
                sumw += weight;
            }
            if (sumw) sum /= sumw;

            item->find<Log>()->info("   Metallicity: " + StringUtils::toString(sum));
        }
    }

    void probeMedium(const SimulationItem* item, int index, const ImportedMedium* medium)
    {
        item->find<Log>()->info("== " + medium->type() + " " + std::to_string(index) + " ==");
        auto snapshot = medium->snapshot();

        EntityCollection entities;
        snapshot->getEntities(entities, Position(), Direction(1,0,0));

        if (snapshot->hasTemperature())
        {
            double sum = 0.;
            double sumw = 0.;
            for (const auto& entity : entities)
            {
                double value = snapshot->temperature(entity.first);
                double weight = snapshot->density(entity.first) * entity.second;
                sum += value * weight;
                sumw += weight;
            }
            if (sumw) sum /= sumw;

            item->find<Log>()->info("   Temperature: " + StringUtils::toString(sum) + " K");
        }
    }
}

////////////////////////////////////////////////////////////////////

void InputModelProbe::probeSetup()
{
    // --- sources ---
    const auto& sources = find<SourceSystem>()->sources();
    int numSources = sources.size();
    for (int i = 0; i != numSources; ++i)
    {
        auto source = dynamic_cast<ImportedSource*>(sources[i]);
        if (source) probeSource(this, i, source);
    }

    // --- media ---
    if (find<Configuration>()->hasMedium())
    {
        const auto& media = find<MediumSystem>()->media();
        int numMedia = media.size();
        for (int i = 0; i != numMedia; ++i)
        {
            auto medium = dynamic_cast<ImportedMedium*>(media[i]);
            if (medium) probeMedium(this, i, medium);
        }
    }
}

////////////////////////////////////////////////////////////////////
