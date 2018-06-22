/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ParticleSnapshot.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "SmoothingKernel.hpp"
#include "TextInFile.hpp"

////////////////////////////////////////////////////////////////////

void ParticleSnapshot::open(const SimulationItem* item, string filename, string description)
{
    Snapshot::open(item, filename, description);
    importPosition();
    importSize();
}

////////////////////////////////////////////////////////////////////

void ParticleSnapshot::readAndClose()
{
    // read the particle info into memory
    // if the user configured a temperature cutoff, we need to skip the "hot" particles
    int numIgnored = 0;
    if (!hasTemperatureCutoff()) _propv = infile()->readAllRows();
    else
    {
        Array row;
        while (infile()->readRow(row))
        {
            if (row[temperatureIndex()] > maxTemperature()) numIgnored++;
            else _propv.push_back(row);
        }
    }

    // close the file
    Snapshot::readAndClose();

    // log the number of particles
    if (!numIgnored)
    {
        log()->info("  Number of particles: " + std::to_string(_propv.size()));
    }
    else
    {
        log()->info("  Number of particles ignored: " + std::to_string(numIgnored));
        log()->info("  Number of particles retained: " + std::to_string(_propv.size()));
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

    // otherwise find the spatial range of the particles assuming a finite support kernel
    double xmin = + std::numeric_limits<double>::infinity();
    double xmax = - std::numeric_limits<double>::infinity();
    double ymin = + std::numeric_limits<double>::infinity();
    double ymax = - std::numeric_limits<double>::infinity();
    double zmin = + std::numeric_limits<double>::infinity();
    double zmax = - std::numeric_limits<double>::infinity();
    for (const Array& prop : _propv)
    {
        xmin = min(xmin, prop[positionIndex()+0] - prop[sizeIndex()]);
        xmax = max(xmax, prop[positionIndex()+0] + prop[sizeIndex()]);
        ymin = min(ymin, prop[positionIndex()+1] - prop[sizeIndex()]);
        ymax = max(ymax, prop[positionIndex()+1] + prop[sizeIndex()]);
        zmin = min(zmin, prop[positionIndex()+2] - prop[sizeIndex()]);
        zmax = max(zmax, prop[positionIndex()+2] + prop[sizeIndex()]);
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
    return Position(_propv[m][positionIndex()+0], _propv[m][positionIndex()+1], _propv[m][positionIndex()+2]);
}

////////////////////////////////////////////////////////////////////

Vec ParticleSnapshot::velocity(int m) const
{
    return Vec(_propv[m][velocityIndex()+0], _propv[m][velocityIndex()+1], _propv[m][velocityIndex()+2]);
}

////////////////////////////////////////////////////////////////////

double ParticleSnapshot::density(Position bfr) const
{

}

////////////////////////////////////////////////////////////////////

double ParticleSnapshot::mass() const
{

}

////////////////////////////////////////////////////////////////////

Position ParticleSnapshot::generatePosition(int m) const
{
    // get center position and size for this particle
    Position ctr(_propv[m][positionIndex()+0], _propv[m][positionIndex()+1], _propv[m][positionIndex()+2]);
    double h = _propv[m][sizeIndex()];

    // sample random position inside the smoothed unit volume
    double u = _kernel->generateRadius();
    Direction k = random()->direction();

    return Position(ctr + k*u*h);
}

////////////////////////////////////////////////////////////////////

Position ParticleSnapshot::generatePosition() const
{
    int m = NR::locateClip(_cumrhov, random()->uniform());
    return generatePosition(m);
}

////////////////////////////////////////////////////////////////////
