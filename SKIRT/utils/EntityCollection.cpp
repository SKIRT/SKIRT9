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

double EntityCollection::average(std::function<double(int m)> value, std::function<double(int m)> weight)
{
    if (_entities.size() > 1)
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

Vec EntityCollection::average(std::function<Vec(int)> value, std::function<double(int)> weight)
{
    if (_entities.size() > 1)
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
        return sumw > 0. ? sumvw / sumw : Vec();
    }
    else if (_entities.size() == 1)
    {
        return value(_entities.cbegin()->first);
    }
    else
    {
        return Vec();
    }
}

////////////////////////////////////////////////////////////////////
