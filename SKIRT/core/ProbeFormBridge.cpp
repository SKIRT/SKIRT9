/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ProbeFormBridge.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "Form.hpp"
#include "MediumSystem.hpp"
#include "PathSegmentGenerator.hpp"
#include "Probe.hpp"
#include "SpatialGrid.hpp"
#include "SpatialGridPath.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

ProbeFormBridge::ProbeFormBridge(const Probe* probe, const Form* form)
{
    // remember the probe and form being bridged
    _probe = probe;
    _form = form;

    // find the simulation's spatial grid, if present
    if (probe->find<Configuration>()->hasMedium()) _grid = probe->find<MediumSystem>()->grid();

    // find the simulation's units system
    _units = probe->find<Units>();
}

////////////////////////////////////////////////////////////////////

void ProbeFormBridge::writeQuantity(string fileid, string quantity, string projectedQuantity, string description,
                                    string projectedDescription, AddColumnDefinitions addColumnDefinitions,
                                    ScalarValueInCell valueInCell)
{
    _type = Type::GridScalarAccumulated;

    _fileid = fileid;
    _unit = _units->unit(quantity);
    _projectedUnit = _units->unit(projectedQuantity);
    _unitFactor = _units->out(quantity, 1.);
    _projectedUnitFactor = _units->out(projectedQuantity, 1.);
    _description = description;
    _projectedDescription = projectedDescription;
    _axis.resize(0);
    _axisUnit = "1";
    _numValues = 1.;

    _addColumnDefinitions = addColumnDefinitions;
    _scalarValueInCell = valueInCell;

    _form->writeQuantity(this);
}

////////////////////////////////////////////////////////////////////

void ProbeFormBridge::writeQuantity(string fileid, string quantity, string description, string projectedDescription,
                                    AddColumnDefinitions addColumnDefinitions, ScalarValueInCell valueInCell,
                                    WeightInCell weightInCell)
{
    _type = Type::GridScalarAveraged;

    _fileid = fileid;
    _unit = _units->unit(quantity);
    _projectedUnit = _unit;
    _unitFactor = _units->out(quantity, 1.);
    _projectedUnitFactor = _unitFactor;
    _description = description;
    _projectedDescription = projectedDescription;
    _axis.resize(0);
    _axisUnit = "1";
    _numValues = 1.;

    _addColumnDefinitions = addColumnDefinitions;
    _scalarValueInCell = valueInCell;
    _weightInCell = weightInCell;

    _form->writeQuantity(this);
}

////////////////////////////////////////////////////////////////////

void ProbeFormBridge::writeQuantity(string fileid, string quantity, string description, string projectedDescription,
                                    AddColumnDefinitions addColumnDefinitions, VectorValueInCell valueInCell,
                                    WeightInCell weightInCell)
{
    _type = Type::GridVectorAveraged;

    _fileid = fileid;
    _unit = _units->unit(quantity);
    _projectedUnit = _unit;
    _unitFactor = _units->out(quantity, 1.);
    _projectedUnitFactor = _unitFactor;
    _description = description;
    _projectedDescription = projectedDescription;
    _axis.resize(3);
    _axis[0] = 1.;
    _axis[1] = 2.;
    _axis[2] = 3.;
    _axisUnit = "1";
    _numValues = 3.;

    _addColumnDefinitions = addColumnDefinitions;
    _vectorValueInCell = valueInCell;
    _weightInCell = weightInCell;

    _form->writeQuantity(this);
}

////////////////////////////////////////////////////////////////////

void ProbeFormBridge::writeQuantity(string fileid, string unit, string description, string projectedDescription,
                                    const Array& axis, string axisUnit, AddColumnDefinitions addColumnDefinitions,
                                    CompoundValueInCell valueInCell, WeightInCell weightInCell)
{
    _type = Type::GridCompoundAveraged;

    _fileid = fileid;
    _unit = unit;
    _projectedUnit = _unit;
    _unitFactor = 1.;
    _projectedUnitFactor = 1.;
    _description = description;
    _projectedDescription = projectedDescription;
    _axis = axis;
    _axisUnit = axisUnit;
    _numValues = axis.size();

    _addColumnDefinitions = addColumnDefinitions;
    _compoundValueInCell = valueInCell;
    _weightInCell = weightInCell;

    _form->writeQuantity(this);
}

////////////////////////////////////////////////////////////////////

void ProbeFormBridge::writeQuantity(string fileid, string quantity, string projectedQuantity, string description,
                                    string projectedDescription, AddColumnDefinitions addColumnDefinitions,
                                    ScalarValueAtPosition valueAtPosition, ScalarValueAlongPath valueAlongPath)
{
    _type = Type::InputScalar;

    _fileid = fileid;
    _unit = _units->unit(quantity);
    _projectedUnit = _units->unit(projectedQuantity);
    _unitFactor = _units->out(quantity, 1.);
    _projectedUnitFactor = _units->out(projectedQuantity, 1.);
    _description = description;
    _projectedDescription = projectedDescription;
    _axis.resize(0);
    _axisUnit = "1";
    _numValues = 1.;

    _addColumnDefinitions = addColumnDefinitions;
    _scalarValueAtPosition = valueAtPosition;
    _scalarValueAlongPath = valueAlongPath;

    _form->writeQuantity(this);
}

////////////////////////////////////////////////////////////////////

void ProbeFormBridge::writeQuantity(string fileid, string quantity, string projectedQuantity, string description,
                                    string projectedDescription, AddColumnDefinitions addColumnDefinitions,
                                    VectorValueAtPosition valueAtPosition, VectorValueAlongPath valueAlongPath)
{
    _type = Type::InputVector;

    _fileid = fileid;
    _unit = _units->unit(quantity);
    _projectedUnit = _units->unit(projectedQuantity);
    _unitFactor = _units->out(quantity, 1.);
    _projectedUnitFactor = _units->out(projectedQuantity, 1.);
    _description = description;
    _projectedDescription = projectedDescription;
    _axis.resize(3);
    _axis[0] = 1.;
    _axis[1] = 2.;
    _axis[2] = 3.;
    _axisUnit = "1";
    _numValues = 3.;

    _addColumnDefinitions = addColumnDefinitions;
    _vectorValueAtPosition = valueAtPosition;
    _vectorValueAlongPath = valueAlongPath;

    _form->writeQuantity(this);
}

////////////////////////////////////////////////////////////////////

void ProbeFormBridge::writeQuantity(string fileid, string unit, string projectedUnit, string description,
                                    string projectedDescription, const Array& axis, string axisUnit,
                                    AddColumnDefinitions addColumnDefinitions, CompoundValueAtPosition valueAtPosition,
                                    CompoundValueAlongPath valueAlongPath)
{
    _type = Type::InputCompound;

    _fileid = fileid;
    _unit = unit;
    _projectedUnit = projectedUnit;
    _unitFactor = 1.;
    _projectedUnitFactor = 1.;
    _description = description;
    _projectedDescription = projectedDescription;
    _axis = axis;
    _axisUnit = axisUnit;
    _numValues = 3.;

    _addColumnDefinitions = addColumnDefinitions;
    _compoundValueAtPosition = valueAtPosition;
    _compoundValueAlongPath = valueAlongPath;

    _form->writeQuantity(this);
}

////////////////////////////////////////////////////////////////////

const SimulationItem* ProbeFormBridge::probe() const
{
    return _probe;
}

////////////////////////////////////////////////////////////////////

const SpatialGrid* ProbeFormBridge::grid() const
{
    return _grid;
}

////////////////////////////////////////////////////////////////////

const Units* ProbeFormBridge::units() const
{
    return _units;
}

////////////////////////////////////////////////////////////////////

string ProbeFormBridge::prefix() const
{
    string result = _probe->itemName();
    if (!_fileid.empty()) result += "_" + _fileid;
    return result;
}

////////////////////////////////////////////////////////////////////

string ProbeFormBridge::unit() const
{
    return _unit;
}

////////////////////////////////////////////////////////////////////

string ProbeFormBridge::projectedUnit() const
{
    return _projectedUnit;
}

////////////////////////////////////////////////////////////////////

string ProbeFormBridge::description() const
{
    return _description;
}

////////////////////////////////////////////////////////////////////

string ProbeFormBridge::projectedDescription() const
{
    return _projectedDescription;
}

////////////////////////////////////////////////////////////////////

Array ProbeFormBridge::axis() const
{
    return _axis;
}

////////////////////////////////////////////////////////////////////

string ProbeFormBridge::axisUnit() const
{
    return _axisUnit;
}

////////////////////////////////////////////////////////////////////

int ProbeFormBridge::numValues() const
{
    return _numValues;
}

////////////////////////////////////////////////////////////////////

bool ProbeFormBridge::isVector() const
{
    return _type == Type::GridVectorAveraged || _type == Type::InputVector;
}

////////////////////////////////////////////////////////////////////

void ProbeFormBridge::addColumnDefinitions(TextOutFile& outfile) const
{
    _addColumnDefinitions(outfile);
}

////////////////////////////////////////////////////////////////////

void ProbeFormBridge::valuesInCell(int m, Array& values) const
{
    switch (_type)
    {
        case Type::GridScalarAccumulated:
        case Type::GridScalarAveraged:
        {
            values[0] = _scalarValueInCell(m) * _unitFactor;
            break;
        }
        case Type::GridVectorAveraged:
        {
            Vec v = _vectorValueInCell(m) * _unitFactor;
            values[0] = v.x();
            values[1] = v.y();
            values[2] = v.z();
            break;
        }
        case Type::GridCompoundAveraged:
        {
            values = _compoundValueInCell(m);
            break;
        }
        case Type::InputScalar:
        case Type::InputVector:
        case Type::InputCompound:
        {
            throw FATALERROR("Cannot retrieve values in given spatial cell from input model probe");
        }
    }
}

////////////////////////////////////////////////////////////////////

void ProbeFormBridge::valuesAtPosition(Position bfr, Array& values) const
{
    switch (_type)
    {
        case Type::GridScalarAccumulated:
        case Type::GridScalarAveraged:
        {
            int m = _grid ? _grid->cellIndex(bfr) : -1;
            if (m >= 0)
                values[0] = _scalarValueInCell(m) * _unitFactor;
            else
                values[0] = 0.;
            break;
        }
        case Type::GridVectorAveraged:
        {
            int m = _grid ? _grid->cellIndex(bfr) : -1;
            if (m >= 0)
            {
                Vec v = _vectorValueInCell(m) * _unitFactor;
                values[0] = v.x();
                values[1] = v.y();
                values[2] = v.z();
            }
            else
            {
                values[0] = 0.;
                values[1] = 0.;
                values[2] = 0.;
            }
            break;
        }
        case Type::GridCompoundAveraged:
        {
            int m = _grid ? _grid->cellIndex(bfr) : -1;
            if (m >= 0)
                values = _compoundValueInCell(m);
            else
                values = 0.;
            break;
        }
        case Type::InputScalar:
        {
            values[0] = _scalarValueAtPosition(bfr) * _unitFactor;
            break;
        }
        case Type::InputVector:
        {
            Vec v = _vectorValueAtPosition(bfr) * _unitFactor;
            values[0] = v.x();
            values[1] = v.y();
            values[2] = v.z();
            break;
        }
        case Type::InputCompound:
        {
            values = _compoundValueAtPosition(bfr);
            break;
        }
    }
}

////////////////////////////////////////////////////////////////////

void ProbeFormBridge::valuesAlongPath(Position bfr, Direction bfk, Array& values) const
{
    switch (_type)
    {
        case Type::GridScalarAccumulated:
        {
            values[0] = 0.;
            if (_grid)
            {
                // get a segment generator and initialize the path
                auto generator = _grid->createPathSegmentGenerator();
                SpatialGridPath path(bfr, bfk);
                generator->start(&path);

                // accumulate values along the path
                while (generator->next())
                {
                    int m = generator->m();
                    if (m >= 0) values[0] += generator->ds() * _scalarValueInCell(m);
                }
                values[0] *= _projectedUnitFactor;
            }
            break;
        }
        case Type::GridScalarAveraged:
        {
            values[0] = 0.;
            if (_grid)
            {
                // get a segment generator and initialize the path
                auto generator = _grid->createPathSegmentGenerator();
                SpatialGridPath path(bfr, bfk);
                generator->start(&path);

                // average values along the path using weights provided by probe
                double totalWeight = 0.;
                while (generator->next())
                {
                    int m = generator->m();
                    if (m >= 0)
                    {
                        double weight = _weightInCell(m) * generator->ds();
                        totalWeight += weight;
                        values[0] += weight * _scalarValueInCell(m);
                    }
                }
                values[0] *= _projectedUnitFactor / totalWeight;
            }
            break;
        }
        case Type::GridVectorAveraged:
        {
            Vec v;
            if (_grid)
            {
                // get a segment generator and initialize the path
                auto generator = _grid->createPathSegmentGenerator();
                SpatialGridPath path(bfr, bfk);
                generator->start(&path);

                // average values along the path using weights provided by probe
                double totalWeight = 0.;
                while (generator->next())
                {
                    int m = generator->m();
                    if (m >= 0)
                    {
                        double weight = _weightInCell(m) * generator->ds();
                        totalWeight += weight;
                        v += weight * _vectorValueInCell(m);
                    }
                }
                v *= _projectedUnitFactor / totalWeight;
            }
            values[0] = v.x();
            values[1] = v.y();
            values[2] = v.z();
            break;
        }
        case Type::GridCompoundAveraged:
        {
            values = 0.;
            if (_grid)
            {
                // get a segment generator and initialize the path
                auto generator = _grid->createPathSegmentGenerator();
                SpatialGridPath path(bfr, bfk);
                generator->start(&path);

                // average values along the path using weights provided by probe
                double totalWeight = 0.;
                while (generator->next())
                {
                    int m = generator->m();
                    if (m >= 0)
                    {
                        double weight = _weightInCell(m) * generator->ds();
                        totalWeight += weight;
                        values += weight * _compoundValueInCell(m);
                    }
                }
                values /= totalWeight;
            }
            break;
        }
        case Type::InputScalar:
        {
            values[0] = _scalarValueAlongPath(bfr, bfk) * _projectedUnitFactor;
            break;
        }
        case Type::InputVector:
        {
            Vec v = _vectorValueAlongPath(bfr, bfk) * _projectedUnitFactor;
            values[0] = v.x();
            values[1] = v.y();
            values[2] = v.z();
            break;
        }
        case Type::InputCompound:
        {
            values = _compoundValueAlongPath(bfr, bfk);
            break;
        }
    }
}

////////////////////////////////////////////////////////////////////
