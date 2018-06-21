/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Snapshot.hpp"
#include "TextInFile.hpp"

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

void Snapshot::open(const SimulationItem* item, std::string filename, std::string description)
{
    _infile = new TextInFile(item, filename, description);
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
}

////////////////////////////////////////////////////////////////////
