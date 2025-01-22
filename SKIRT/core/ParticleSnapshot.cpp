/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ParticleSnapshot.hpp"
#include "EntityCollection.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "SmoothingKernel.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"

////////////////////////////////////////////////////////////////////

/** Particle is a simple helper class for working with smoothed particles, usually
    imported from a smoothed particle hydrodynamical (SPH) simulation. A Particle object
    holds the particle's center position, smoothing length, and effective mass (after any
    multipliers have been applied). Isolating these properties in a simple object allows efficient
    storage and retrieval, for example when organizing smoothed particles in a search structure. */
class ParticleSnapshot::Particle
{
public:
    /** The constructor arguments specify the particle attributes: particle index, coordinates of
        the center, smoothing length, and effective mass. */
    Particle(double x, double y, double z, double h, double M) : _x(x), _y(y), _z(z), _h(h), _M(M) {}

    /** This function returns the coordinates of the center of the particle. */
    Vec center() const { return Vec(_x, _y, _z); }

    /** This function returns the smoothing length of the particle. */
    double radius() const { return _h; }

    /** This function returns the bounding box of the particle. */
    Box bounds() const { return Box(_x - _h, _y - _h, _z - _h, _x + _h, _y + _h, _z + _h); }

    /** This function returns the mass of the particle. */
    double mass() const { return _M; }

    /** This function returns the effective volume of the particle. */
    double volume() const { return _h * _h * _h; }

    /** This function returns the effective density of the particle. */
    double density() const { return _M / (_h * _h * _h); }

    /** This function returns the impact radius of the path specified by a starting position and
        direction relative to the particle center (i.e. the minimum distance between the path and
        the center). If the path's starting position is located beyond the impact point (the point
        of minimum distance), the function returns infinity to indicate that the path does not
        intersect the particle. */
    double impact(Vec r, Vec k) const
    {
        // assume that k is normalized so we don't need to divide by its norm
        double s = Vec::dot(k, center() - r);
        if (s < 0.) return std::numeric_limits<double>::infinity();
        return (center() - (r + s * k)).norm();
    }

private:
    // data members received as constructor arguments
    double _x, _y, _z;  // center coordinates
    double _h;          // smoothing length
    double _M;          // total mass
};

////////////////////////////////////////////////////////////////////

ParticleSnapshot::ParticleSnapshot() {}

////////////////////////////////////////////////////////////////////

ParticleSnapshot::~ParticleSnapshot() {}

////////////////////////////////////////////////////////////////////

void ParticleSnapshot::readAndClose()
{
    // read the particle info into memory
    // if the user configured a temperature cutoff, we skip high-temperature particles
    // if the user configured a mass-density policy, we skip zero-mass particles
    int numTempIgnored = 0;
    int numMassIgnored = 0;
    int numBiasIgnored = 0;
    Array row;
    while (infile()->readRow(row))
    {
        if (useTemperatureCutoff() && row[temperatureIndex()] > maxTemperature())
            numTempIgnored++;
        else if (hasMassDensityPolicy() && row[massIndex()] == 0.)
            numMassIgnored++;
        else if (hasBias() && row[biasIndex()] == 0.)
            numBiasIgnored++;
        else
            _propv.push_back(row);
    }

    // close the file
    close();

    // log the number of particles
    if (!numTempIgnored && !numMassIgnored && !numBiasIgnored)
    {
        log()->info("  Number of particles: " + std::to_string(_propv.size()));
    }
    else
    {
        if (numTempIgnored)
            log()->info("  Number of high-temperature particles ignored: " + std::to_string(numTempIgnored));
        if (numMassIgnored) log()->info("  Number of zero-mass particles ignored: " + std::to_string(numMassIgnored));
        if (numBiasIgnored) log()->info("  Number of zero-bias particles ignored: " + std::to_string(numBiasIgnored));
        log()->info("  Number of particles retained: " + std::to_string(_propv.size()));
    }

    // if a mass density policy has been set, calculate masses and densities for all cells
    if (hasMassDensityPolicy())
    {
        // build a list of compact smoothed particle objects that we can organize in a search structure
        double totalOriginalMass = 0;
        double totalMetallicMass = 0;
        double totalEffectiveMass = 0;
        int numParticles = _propv.size();
        _pv.reserve(numParticles);
        for (int m = 0; m != numParticles; ++m)
        {
            const Array& prop = _propv[m];

            double originalMass = prop[massIndex()];
            double metallicMass = originalMass * (useMetallicity() ? prop[metallicityIndex()] : 1.);
            double effectiveMass = metallicMass * multiplier();

            _pv.emplace_back(prop[positionIndex() + 0], prop[positionIndex() + 1], prop[positionIndex() + 2],
                             prop[sizeIndex()], effectiveMass);

            totalOriginalMass += originalMass;
            totalMetallicMass += metallicMass;
            totalEffectiveMass += effectiveMass;
        }

        // log mass statistics
        logMassStatistics(0, totalOriginalMass, totalMetallicMass, totalEffectiveMass);

        // if one of the total masses is negative, suppress the complete mass distribution
        if (totalOriginalMass < 0 || totalMetallicMass < 0 || totalEffectiveMass < 0)
        {
            log()->warning("  Total imported mass is negative; suppressing the complete mass distribution");
            _propv.clear();
            _pv.clear();
            totalEffectiveMass = 0;
        }

        // remember the effective mass
        _mass = totalEffectiveMass;

        // construct a vector with the normalized cumulative particle densities
        if (!_pv.empty()) NR::cdf(_cumrhov, _pv.size(), [this](int i) { return _pv[i].mass(); });
    }

    // if needed, build a list of compact particles without masses that we can organize in a search structure
    if (!hasMassDensityPolicy() && needGetEntities())
    {
        int numParticles = _propv.size();
        _pv.reserve(numParticles);
        for (int m = 0; m != numParticles; ++m)
        {
            const Array& prop = _propv[m];
            _pv.emplace_back(prop[positionIndex() + 0], prop[positionIndex() + 1], prop[positionIndex() + 2],
                             prop[sizeIndex()], 0.);
        }
    }

    // if needed, construct a search structure for the particles
    if (hasMassDensityPolicy() || needGetEntities())
    {
        log()->info("Constructing search grid for " + std::to_string(_pv.size()) + " particles...");
        auto bounds = [this](int m) { return _pv[m].bounds(); };
        auto intersects = [this](int m, const Box& box) { return box.intersects(_pv[m].center(), _pv[m].radius()); };
        _search.loadEntities(_pv.size(), bounds, intersects);

        int nb = _search.numBlocks();
        log()->info("  Number of blocks in grid: " + std::to_string(nb * nb * nb) + " (" + std::to_string(nb) + "^3)");
        log()->info("  Smallest number of particles per block: " + std::to_string(_search.minEntitiesPerBlock()));
        log()->info("  Largest  number of particles per block: " + std::to_string(_search.maxEntitiesPerBlock()));
        log()->info("  Average  number of particles per block: "
                    + StringUtils::toString(_search.avgEntitiesPerBlock(), 'f', 1));
    }
}

////////////////////////////////////////////////////////////////////

void ParticleSnapshot::setSmoothingKernel(const SmoothingKernel* kernel)
{
    _kernel = kernel;
}

////////////////////////////////////////////////////////////////////

Box ParticleSnapshot::extent() const
{
    // if there are no particles, return an empty box
    if (_propv.empty()) return Box();

    // if there is a search structure, ask it to return the extent (it is already calculated)
    if (_search.numBlocks()) return _search.extent();

    // otherwise find the spatial range of the particles assuming a finite support kernel
    double xmin = +std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    double ymin = +std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();
    double zmin = +std::numeric_limits<double>::infinity();
    double zmax = -std::numeric_limits<double>::infinity();
    for (const Array& prop : _propv)
    {
        xmin = min(xmin, prop[positionIndex() + 0] - prop[sizeIndex()]);
        xmax = max(xmax, prop[positionIndex() + 0] + prop[sizeIndex()]);
        ymin = min(ymin, prop[positionIndex() + 1] - prop[sizeIndex()]);
        ymax = max(ymax, prop[positionIndex() + 1] + prop[sizeIndex()]);
        zmin = min(zmin, prop[positionIndex() + 2] - prop[sizeIndex()]);
        zmax = max(zmax, prop[positionIndex() + 2] + prop[sizeIndex()]);
    }
    return Box(xmin, ymin, zmin, xmax, ymax, zmax);
}

////////////////////////////////////////////////////////////////////

int ParticleSnapshot::numEntities() const
{
    return _propv.size();
}

////////////////////////////////////////////////////////////////////

double ParticleSnapshot::volume(int m) const
{
    return _pv[m].volume();
}

////////////////////////////////////////////////////////////////////

double ParticleSnapshot::density(int m) const
{
    return _pv[m].density();
}

////////////////////////////////////////////////////////////////////

double ParticleSnapshot::density(Position bfr) const
{
    double sum = 0.;
    for (int m : _search.entitiesFor(bfr))
    {
        double u = (bfr - _pv[m].center()).norm() / _pv[m].radius();
        sum += _kernel->density(u) * _pv[m].density();
    }
    return sum > 0. ? sum : 0.;  // guard against negative densities
}

////////////////////////////////////////////////////////////////////

double ParticleSnapshot::mass() const
{
    return _mass;
}

////////////////////////////////////////////////////////////////////

Position ParticleSnapshot::position(int m) const
{
    return Position(_propv[m][positionIndex() + 0], _propv[m][positionIndex() + 1], _propv[m][positionIndex() + 2]);
}

////////////////////////////////////////////////////////////////////

Position ParticleSnapshot::generatePosition(int m) const
{
    // get center position and size for this particle
    Position rc(_propv[m][positionIndex() + 0], _propv[m][positionIndex() + 1], _propv[m][positionIndex() + 2]);
    double h = _propv[m][sizeIndex()];

    // sample random position inside the smoothed unit volume
    double u = _kernel->generateRadius();
    Direction k = random()->direction();

    return Position(rc + k * u * h);
}

////////////////////////////////////////////////////////////////////

Position ParticleSnapshot::generatePosition() const
{
    // if there are no particles, return the origin
    if (_propv.empty()) return Position();

    // select a particle according to its mass contribution
    int m = NR::locateClip(_cumrhov, random()->uniform());

    return generatePosition(m);
}

////////////////////////////////////////////////////////////////////

const Array& ParticleSnapshot::properties(int m) const
{
    return _propv[m];
}

////////////////////////////////////////////////////////////////////

void ParticleSnapshot::getEntities(EntityCollection& entities, Position bfr) const
{
    entities.clear();
    for (int m : _search.entitiesFor(bfr))
    {
        double u = (bfr - _pv[m].center()).norm() / _pv[m].radius();
        entities.add(m, _kernel->density(u));
    };
}

////////////////////////////////////////////////////////////////////

void ParticleSnapshot::getEntities(EntityCollection& entities, Position bfr, Direction bfk) const
{
    entities.clear();
    for (int m : _search.entitiesFor(bfr, bfk))
    {
        double h = _pv[m].radius();
        double q = _pv[m].impact(bfr, bfk) / h;
        entities.add(m, _kernel->columnDensity(q) * h);
    };
}

////////////////////////////////////////////////////////////////////
