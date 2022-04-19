/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "EntityCollection.hpp"

////////////////////////////////////////////////////////////////////

EntityCollection::EntityCollection() {}

////////////////////////////////////////////////////////////////////

void EntityCollection::clear()
{
    _entities.clear();
}

////////////////////////////////////////////////////////////////////

void EntityCollection::add(int m, double w)
{
    if (m >= 0 && w > 0. && std::isfinite(w)) _entities.emplace(m, w);
}

////////////////////////////////////////////////////////////////////

void EntityCollection::addSingle(int m)
{
    _entities.clear();
    if (m >= 0) _entities.emplace(m, 1.);
}

////////////////////////////////////////////////////////////////////

double EntityCollection::average(std::function<double(int)> value)
{
    if (_entities.size() > 1)
    {
        double sumvw = 0.;
        double sumw = 0.;
        for (const auto& entity : _entities)
        {
            sumvw += value(entity.first) * entity.second;
            sumw += entity.second;
        }
        return sumw > 0. ? sumvw / sumw : 0.;
    }
    else if (_entities.size() == 1)
    {
        return value(_entities.cbegin()->first);
    }
    else
    {
        return 0.;
    }
}

////////////////////////////////////////////////////////////////////

double EntityCollection::accumulate(std::function<double(int)> value)
{
    double sumvw = 0.;
    for (const auto& entity : _entities)
    {
        sumvw += value(entity.first) * entity.second;
    }
    return sumvw;
}

////////////////////////////////////////////////////////////////////
