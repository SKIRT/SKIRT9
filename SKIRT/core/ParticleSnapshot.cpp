/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ParticleSnapshot.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "SmoothedParticleGrid.hpp"
#include "SmoothingKernel.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

ParticleSnapshot::~ParticleSnapshot()
{
    delete _grid;
}

////////////////////////////////////////////////////////////////////

void ParticleSnapshot::readAndClose()
{
    // read the particle info into memory
    // if the user configured a temperature cutoff, we skip high-temperature particles
    // if the user configured a mass-density policy, we skip zero-mass particles
    int numTempIgnored = 0;
    int numMassIgnored = 0;
    Array row;
    while (infile()->readRow(row))
    {
        if (useTemperatureCutoff() && row[temperatureIndex()] > maxTemperature())
            numTempIgnored++;
        else if (hasMassDensityPolicy() && row[massIndex()] == 0)
            numMassIgnored++;
        else
            _propv.push_back(row);
    }

    // close the file
    Snapshot::readAndClose();

    // log the number of particles
    if (!numTempIgnored && !numMassIgnored)
    {
        log()->info("  Number of particles: " + std::to_string(_propv.size()));
    }
    else
    {
        if (numTempIgnored)
            log()->info("  Number of high-temperature particles ignored: " + std::to_string(numTempIgnored));
        if (numMassIgnored) log()->info("  Number of zero-mass particles ignored: " + std::to_string(numMassIgnored));
        log()->info("  Number of particles retained: " + std::to_string(_propv.size()));
    }

    // we can calculate mass and densities only if a policy has been set
    if (!hasMassDensityPolicy()) return;

    // build a list of compact smoothed particle objects that we can organize in a grid
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

        _pv.emplace_back(m, prop[positionIndex() + 0], prop[positionIndex() + 1], prop[positionIndex() + 2],
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
        return;  // abort
    }

    // remember the effective mass
    _mass = totalEffectiveMass;

    // if there are no particles, do not build the special structures for optimizing operations
    if (_pv.empty()) return;

    // construct a 3D-grid over the particle space, and create a list of particles that overlap each grid cell
    int gridsize = max(20, static_cast<int>(pow(_pv.size(), 1. / 3.) / 5));
    string size = std::to_string(gridsize);
    log()->info("Constructing intermediate " + size + "x" + size + "x" + size + " grid for particles...");
    _grid = new SmoothedParticleGrid(_pv, gridsize);
    log()->info("  Smallest number of particles per cell: " + std::to_string(_grid->minParticlesPerCell()));
    log()->info("  Largest  number of particles per cell: " + std::to_string(_grid->maxParticlesPerCell()));
    log()->info("  Average  number of particles per cell: "
                + StringUtils::toString(_grid->totalParticles() / double(gridsize * gridsize * gridsize), 'f', 1));

    // construct a vector with the normalized cumulative particle densities
    NR::cdf(_cumrhov, _pv.size(), [this](int i) { return _pv[i].mass(); });
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

    // if there is a particle grid, ask it to return the extent (it is already calculated)
    if (_grid) return _grid->extent();

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

Position ParticleSnapshot::position(int m) const
{
    return Position(_propv[m][positionIndex() + 0], _propv[m][positionIndex() + 1], _propv[m][positionIndex() + 2]);
}

////////////////////////////////////////////////////////////////////

double ParticleSnapshot::metallicity(int m) const
{
    return _propv[m][metallicityIndex()];
}

////////////////////////////////////////////////////////////////////

double ParticleSnapshot::metallicity(Position bfr) const
{
    const SmoothedParticle* nearestParticle = _grid ? _grid->nearestParticle(bfr) : nullptr;
    return nearestParticle ? metallicity(nearestParticle->index()) : 0.;
}

////////////////////////////////////////////////////////////////////

double ParticleSnapshot::temperature(int m) const
{
    return _propv[m][temperatureIndex()];
}

////////////////////////////////////////////////////////////////////

double ParticleSnapshot::temperature(Position bfr) const
{
    const SmoothedParticle* nearestParticle = _grid ? _grid->nearestParticle(bfr) : nullptr;
    return nearestParticle ? temperature(nearestParticle->index()) : 0.;
}

////////////////////////////////////////////////////////////////////

Vec ParticleSnapshot::velocity(int m) const
{
    return Vec(_propv[m][velocityIndex() + 0], _propv[m][velocityIndex() + 1], _propv[m][velocityIndex() + 2]);
}

////////////////////////////////////////////////////////////////////

Vec ParticleSnapshot::velocity(Position bfr) const
{
    const SmoothedParticle* nearestParticle = _grid ? _grid->nearestParticle(bfr) : nullptr;
    return nearestParticle ? velocity(nearestParticle->index()) : Vec();
}

////////////////////////////////////////////////////////////////////

double ParticleSnapshot::velocityDispersion(int m) const
{
    return _propv[m][velocityDispersionIndex()];
}

////////////////////////////////////////////////////////////////////

double ParticleSnapshot::velocityDispersion(Position bfr) const
{
    const SmoothedParticle* nearestParticle = _grid ? _grid->nearestParticle(bfr) : nullptr;
    return nearestParticle ? velocityDispersion(nearestParticle->index()) : 0.;
}

////////////////////////////////////////////////////////////////////

Vec ParticleSnapshot::magneticField(int m) const
{
    return Vec(_propv[m][magneticFieldIndex() + 0], _propv[m][magneticFieldIndex() + 1],
               _propv[m][magneticFieldIndex() + 2]);
}

////////////////////////////////////////////////////////////////////

Vec ParticleSnapshot::magneticField(Position bfr) const
{
    const SmoothedParticle* nearestParticle = _grid ? _grid->nearestParticle(bfr) : nullptr;
    return nearestParticle ? magneticField(nearestParticle->index()) : Vec();
}

////////////////////////////////////////////////////////////////////

void ParticleSnapshot::parameters(int m, Array& params) const
{
    int n = numParameters();
    params.resize(n);
    for (int i = 0; i != n; ++i) params[i] = _propv[m][parametersIndex() + i];
}

////////////////////////////////////////////////////////////////////

void ParticleSnapshot::parameters(Position bfr, Array& params) const
{
    const SmoothedParticle* nearestParticle = _grid ? _grid->nearestParticle(bfr) : nullptr;
    if (nearestParticle)
        parameters(nearestParticle->index(), params);
    else
        params.resize(numParameters());
}

////////////////////////////////////////////////////////////////////

double ParticleSnapshot::density(Position bfr) const
{
    double sum = 0.;
    if (_grid)
        for (const SmoothedParticle* p : _grid->particlesFor(bfr))
        {
            double h = p->radius();
            double u = (bfr - p->center()).norm() / h;
            sum += _kernel->density(u) * p->mass() / (h * h * h);
        }
    return sum > 0. ? sum : 0.;  // guard against negative densities
}

////////////////////////////////////////////////////////////////////

double ParticleSnapshot::mass() const
{
    return _mass;
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
