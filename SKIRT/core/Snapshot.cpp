/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Snapshot.hpp"
#include "TextInFile.hpp"
#include "Log.hpp"
#include "Random.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

Snapshot::Snapshot()
{
}

////////////////////////////////////////////////////////////////////

Snapshot::~Snapshot()
{
    delete _infile;
}

////////////////////////////////////////////////////////////////////

void Snapshot::open(const SimulationItem* item, string filename, string description)
{
    _infile = new TextInFile(item, filename, description);
    _log = item->find<Log>();
    _units = item->find<Units>();
    _random = item->find<Random>();
}

////////////////////////////////////////////////////////////////////

void Snapshot::readAndClose()
{
    delete _infile;
    _infile = nullptr;
}

////////////////////////////////////////////////////////////////////

void Snapshot::importPosition()
{
    _positionIndex = _nextIndex++;
    _infile->addColumn("position x", "length", "pc");
    _infile->addColumn("position y", "length", "pc");
    _infile->addColumn("position z", "length", "pc");
}

////////////////////////////////////////////////////////////////////

void Snapshot::importSize()
{
    _sizeIndex = _nextIndex++;
    _infile->addColumn("size h", "length", "pc");
}

////////////////////////////////////////////////////////////////////

void Snapshot::importDensity()
{
    _densityIndex = _nextIndex++;
    _infile->addColumn("mass density", "massvolumedensity", "Msun/pc3");
}

////////////////////////////////////////////////////////////////////

void Snapshot::importMass()
{
    _massIndex = _nextIndex++;
    _infile->addColumn("mass", "mass", "Msun");
}

////////////////////////////////////////////////////////////////////

void Snapshot::importMetallicity()
{
    _metallicityIndex = _nextIndex++;
    _infile->addColumn("metallicity");
}

////////////////////////////////////////////////////////////////////

void Snapshot::importTemperature()
{
    _temperatureIndex = _nextIndex++;
    _infile->addColumn("temperature", "temperature", "K");
}

////////////////////////////////////////////////////////////////////

void Snapshot::importVelocity()
{
    _velocityIndex = _nextIndex++;
    _infile->addColumn("velocity x", "velocity", "km/s");
    _infile->addColumn("velocity y", "velocity", "km/s");
    _infile->addColumn("velocity z", "velocity", "km/s");
}

////////////////////////////////////////////////////////////////////

void Snapshot::importParameters(const vector<SnapshotParameter>& parameters)
{
    _parametersIndex = _nextIndex++;
    _numParameters = parameters.size();
    for (const auto& p : parameters)
        _infile->addColumn(p.description(), p.quantity(), p.defaultUnit());
}

////////////////////////////////////////////////////////////////////

void Snapshot::setMassDensityPolicy(double multiplier, double maxTemperature)
{
    _multiplier = multiplier;
    _maxTemperature = maxTemperature;
    _hasDensityPolicy = true;
}

////////////////////////////////////////////////////////////////////

double Snapshot::volume() const
{
    return extent().volume();
}

////////////////////////////////////////////////////////////////////

namespace
{
    // the number of samples used for integrating the surface densities
    const int NSAMPLES = 10000;
}

////////////////////////////////////////////////////////////////////

double Snapshot::SigmaX() const
{
    // determine a small value relative to the domain extent;
    // we integrate along a small offset from the axis to avoid cell borders
    double eps = 1e-12 * extent().widths().norm();

    double sum = 0;
    double xmin = extent().xmin();
    double xmax = extent().xmax();
    for (int k = 0; k < NSAMPLES; k++)
    {
        sum += density(Position(xmin + k*(xmax-xmin)/NSAMPLES, eps, eps));
    }
    return (sum/NSAMPLES)*(xmax-xmin);
}

////////////////////////////////////////////////////////////////////

double Snapshot::SigmaY() const
{
    // determine a small value relative to the domain extent;
    // we integrate along a small offset from the axis to avoid cell borders
    double eps = 1e-12 * extent().widths().norm();

    double sum = 0;
    double ymin = extent().ymin();
    double ymax = extent().ymax();
    for (int k = 0; k < NSAMPLES; k++)
    {
        sum += density(Position(eps, ymin + k*(ymax-ymin)/NSAMPLES, eps));
    }
    return (sum/NSAMPLES)*(ymax-ymin);
}

////////////////////////////////////////////////////////////////////

double Snapshot::SigmaZ() const
{
    // determine a small value relative to the domain extent;
    // we integrate along a small offset from the axis to avoid cell borders
    double eps = 1e-12 * extent().widths().norm();

    double sum = 0;
    double zmin = extent().zmin();
    double zmax = extent().zmax();
    for (int k = 0; k < NSAMPLES; k++)
    {
        sum += density(Position(eps, eps, zmin + k*(zmax-zmin)/NSAMPLES));
    }
    return (sum/NSAMPLES)*(zmax-zmin);
}

////////////////////////////////////////////////////////////////////
