/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Snapshot.hpp"
#include "EntityCollection.hpp"
#include "Log.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

Snapshot::Snapshot() {}

////////////////////////////////////////////////////////////////////

Snapshot::~Snapshot()
{
    delete _infile;
}

////////////////////////////////////////////////////////////////////

void Snapshot::open(const SimulationItem* item, string filename, string description)
{
    _infile = new TextInFile(item, filename, description);
    setContext(item);
}

////////////////////////////////////////////////////////////////////

void Snapshot::readAndClose()
{
    delete _infile;
    _infile = nullptr;
}

////////////////////////////////////////////////////////////////////

void Snapshot::setContext(const SimulationItem* item)
{
    _log = item->find<Log>();
    _units = item->find<Units>();
    _random = item->find<Random>();
}

////////////////////////////////////////////////////////////////////

void Snapshot::useColumns(string columns)
{
    _infile->useColumns(columns);
}

////////////////////////////////////////////////////////////////////

void Snapshot::importPosition()
{
    _positionIndex = _nextIndex;
    _nextIndex += 3;
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

void Snapshot::importBox()
{
    _boxIndex = _nextIndex;
    _nextIndex += 6;
    _infile->addColumn("box xmin", "length", "pc");
    _infile->addColumn("box ymin", "length", "pc");
    _infile->addColumn("box zmin", "length", "pc");
    _infile->addColumn("box xmax", "length", "pc");
    _infile->addColumn("box ymax", "length", "pc");
    _infile->addColumn("box zmax", "length", "pc");
}

////////////////////////////////////////////////////////////////////

void Snapshot::importMassDensity()
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

void Snapshot::importNumberDensity()
{
    _densityIndex = _nextIndex++;
    _infile->addColumn("number density", "numbervolumedensity", "1/cm3");
    _holdsNumber = true;
}

////////////////////////////////////////////////////////////////////

void Snapshot::importNumber()
{
    _massIndex = _nextIndex++;
    _infile->addColumn("number");
    _holdsNumber = true;
}

////////////////////////////////////////////////////////////////////

void Snapshot::importCurrentMass()
{
    _currentMassIndex = _nextIndex++;
    _infile->addColumn("current mass", "mass", "Msun");
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
    _velocityIndex = _nextIndex;
    _nextIndex += 3;
    _infile->addColumn("velocity x", "velocity", "km/s");
    _infile->addColumn("velocity y", "velocity", "km/s");
    _infile->addColumn("velocity z", "velocity", "km/s");
}

void Snapshot::importVelocityDispersion()
{
    _velocityDispersionIndex = _nextIndex++;
    _infile->addColumn("velocity dispersion", "velocity", "km/s");
}

void Snapshot::importMagneticField()
{
    _magneticFieldIndex = _nextIndex;
    _nextIndex += 3;
    _infile->addColumn("magnetic field x", "magneticfield", "uG");
    _infile->addColumn("magnetic field y", "magneticfield", "uG");
    _infile->addColumn("magnetic field z", "magneticfield", "uG");
}

////////////////////////////////////////////////////////////////////

void Snapshot::importParameters(const vector<SnapshotParameter>& parameters)
{
    _parametersIndex = _nextIndex;
    _numParameters = parameters.size();
    for (const auto& param : parameters)
    {
        _infile->addColumn(param.description(), param.quantity(), param.defaultUnit());

        // discover standard parameter types and set the corresponding index so these values can be probed
        switch (param.identifier())
        {
            case SnapshotParameter::Identifier::InitialMass: _initialMassIndex = _nextIndex; break;
            case SnapshotParameter::Identifier::CurrentMass: _currentMassIndex = _nextIndex; break;
            case SnapshotParameter::Identifier::Metallicity: _metallicityIndex = _nextIndex; break;
            case SnapshotParameter::Identifier::Age: _ageIndex = _nextIndex; break;
            case SnapshotParameter::Identifier::Temperature: _temperatureIndex = _nextIndex; break;
            case SnapshotParameter::Identifier::Custom: break;
        }
        _nextIndex++;
    }
}

////////////////////////////////////////////////////////////////////

void Snapshot::setMassDensityPolicy(double multiplier, double maxTemperature, bool useMetallicity)
{
    _multiplier = multiplier;
    _maxTemperature = maxTemperature;
    _useMetallicity = useMetallicity;
    _hasDensityPolicy = true;
}

////////////////////////////////////////////////////////////////////

void Snapshot::setNeedGetEntities()
{
    _needGetEntities = true;
}

////////////////////////////////////////////////////////////////////

void Snapshot::logMassStatistics(int numIgnored, double totalOriginalMass, double totalMetallicMass,
                                 double totalEffectiveMass)
{
    auto massUnit = holdsNumber() ? "" : units()->umass();
    auto massToString = [this](double mass) {
        return StringUtils::toString(holdsNumber() ? mass : units()->omass(mass), 'e', 4);
    };
    if (numIgnored) log()->info("  Ignored mass in " + std::to_string(numIgnored) + " high-temperature cells");
    log()->info("  Total original mass : " + massToString(totalOriginalMass) + " " + massUnit);
    if (useMetallicity()) log()->info("  Total metallic mass : " + massToString(totalMetallicMass) + " " + massUnit);
    log()->info("  Total effective mass: " + massToString(totalEffectiveMass) + " " + massUnit);
}

////////////////////////////////////////////////////////////////////

double Snapshot::volume() const
{
    return extent().volume();
}

////////////////////////////////////////////////////////////////////

double Snapshot::initialMass(int m) const
{
    return properties(m)[initialMassIndex()];
}

////////////////////////////////////////////////////////////////////

double Snapshot::currentMass(int m) const
{
    return properties(m)[currentMassIndex()];
}

////////////////////////////////////////////////////////////////////

double Snapshot::metallicity(int m) const
{
    return properties(m)[metallicityIndex()];
}

////////////////////////////////////////////////////////////////////

double Snapshot::metallicity(Position bfr) const
{
    thread_local EntityCollection entities;  // can be reused for all queries in a given execution thread
    getEntities(entities, bfr);
    return entities.averageValue([this](int m) { return metallicity(m); }, [this](int m) { return currentMass(m); });
}

////////////////////////////////////////////////////////////////////

double Snapshot::age(int m) const
{
    return properties(m)[ageIndex()];
}

////////////////////////////////////////////////////////////////////

double Snapshot::temperature(int m) const
{
    return properties(m)[temperatureIndex()];
}

////////////////////////////////////////////////////////////////////

double Snapshot::temperature(Position bfr) const
{
    thread_local EntityCollection entities;  // can be reused for all queries in a given execution thread
    getEntities(entities, bfr);
    return entities.averageValue([this](int m) { return temperature(m); }, [this](int m) { return currentMass(m); });
}

////////////////////////////////////////////////////////////////////

Vec Snapshot::velocity(int m) const
{
    const auto& propv = properties(m);
    return Vec(propv[velocityIndex() + 0], propv[velocityIndex() + 1], propv[velocityIndex() + 2]);
}

////////////////////////////////////////////////////////////////////

Vec Snapshot::velocity(Position bfr) const
{
    thread_local EntityCollection entities;  // can be reused for all queries in a given execution thread
    getEntities(entities, bfr);
    return entities.averageValue([this](int m) { return velocity(m); }, [this](int m) { return currentMass(m); });
}

////////////////////////////////////////////////////////////////////

double Snapshot::velocityDispersion(int m) const
{
    return properties(m)[velocityDispersionIndex()];
}

////////////////////////////////////////////////////////////////////

Vec Snapshot::magneticField(int m) const
{
    const auto& propv = properties(m);
    return Vec(propv[magneticFieldIndex() + 0], propv[magneticFieldIndex() + 1], propv[magneticFieldIndex() + 2]);
}

////////////////////////////////////////////////////////////////////

Vec Snapshot::magneticField(Position bfr) const
{
    thread_local EntityCollection entities;  // can be reused for all queries in a given execution thread
    getEntities(entities, bfr);
    return entities.averageValue([this](int m) { return magneticField(m); }, [this](int m) { return currentMass(m); });
}

////////////////////////////////////////////////////////////////////

void Snapshot::parameters(int m, Array& params) const
{
    int n = numParameters();
    params.resize(n);
    const auto& propv = properties(m);
    for (int i = 0; i != n; ++i) params[i] = propv[parametersIndex() + i];
}

////////////////////////////////////////////////////////////////////

void Snapshot::parameters(Position bfr, Array& params) const
{
    thread_local EntityCollection entities;  // can be reused for all queries in a given execution thread
    getEntities(entities, bfr);

    // look for the entity with the highest weight
    double wmax = 0.;
    int m = -1;
    for (const auto& entity : entities)
    {
        if (entity.second > wmax)
        {
            wmax = entity.second;
            m = entity.first;
        }
    }

    // if found, use it
    if (m >= 0)
        parameters(m, params);
    else
        params.resize(numParameters());
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
        sum += density(Position(xmin + k * (xmax - xmin) / NSAMPLES, eps, eps));
    }
    return (sum / NSAMPLES) * (xmax - xmin);
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
        sum += density(Position(eps, ymin + k * (ymax - ymin) / NSAMPLES, eps));
    }
    return (sum / NSAMPLES) * (ymax - ymin);
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
        sum += density(Position(eps, eps, zmin + k * (zmax - zmin) / NSAMPLES));
    }
    return (sum / NSAMPLES) * (zmax - zmin);
}

////////////////////////////////////////////////////////////////////
