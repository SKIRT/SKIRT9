/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "EntityCollection.hpp"
#include "Vec.hpp"

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

std::pair<double, double> EntityCollection::average(std::function<double(int m)> value,
                                                    std::function<double(int m)> weight)
{
    double sumvw = 0.;
    double sumw = 0.;
    for (const auto& entity : _entities)
    {
        double v = value(entity.first);
        double w = weight(entity.first) * entity.second;
        sumvw += v * w;
        sumw += w;
    }
    return std::make_pair(sumvw, sumw);
}

////////////////////////////////////////////////////////////////////

std::pair<Vec, double> EntityCollection::average(std::function<Vec(int m)> value, std::function<double(int m)> weight)
{
    Vec sumvw;
    double sumw = 0.;
    for (const auto& entity : _entities)
    {
        Vec v = value(entity.first);
        double w = weight(entity.first) * entity.second;
        sumvw += v * w;
        sumw += w;
    }
    return std::make_pair(sumvw, sumw);
}

////////////////////////////////////////////////////////////////////
