/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ProbeFormBridge.hpp"
#include "Configuration.hpp"
#include "Direction.hpp"
#include "FatalError.hpp"
#include "Form.hpp"
#include "MediumSystem.hpp"
#include "PathSegmentGenerator.hpp"
#include "Position.hpp"
#include "Probe.hpp"
#include "SpatialGrid.hpp"
#include "SpatialGridPath.hpp"

////////////////////////////////////////////////////////////////////

ProbeFormBridge::ProbeFormBridge(const Probe* probe, const Form* form)
{
    _probe = probe;
    _form = form;

    // find the simulation's spatial grid, if present
    if (_probe->find<Configuration>()->hasMedium()) _grid = _probe->find<MediumSystem>()->grid();
}

////////////////////////////////////////////////////////////////////

void ProbeFormBridge::setQuantity(int numValues, string quantity, string projectedQuantity, string description,
                                  string projectedDescription, AddColumnHeaders addColumnHeaders)
{
    _numValues = numValues;
    _quantity = quantity;
    _projectedQuantity = projectedQuantity;
    _description = description;
    _projectedDescription = projectedDescription;
    _addColumnHeaders = addColumnHeaders;
}

////////////////////////////////////////////////////////////////////

void ProbeFormBridge::writeFile(string fileid, ValuesInCell valuesInCell, WeightInCell weightInCell)
{
    _fileid = fileid;
    _valuesInCell = valuesInCell;
    _weightInCell = weightInCell;
    _valuesAtPosition = nullptr;
    _valuesAlongPath = nullptr;

    _form->writeFile(this);
}

////////////////////////////////////////////////////////////////////

void ProbeFormBridge::writeFile(string fileid, ValuesAtPosition valuesAtPosition, ValuesAlongPath valuesAlongPath)
{
    _fileid = fileid;
    _valuesInCell = nullptr;
    _weightInCell = nullptr;
    _valuesAtPosition = valuesAtPosition;
    _valuesAlongPath = valuesAlongPath;

    _form->writeFile(this);
}

////////////////////////////////////////////////////////////////////

const Probe* ProbeFormBridge::probe() const
{
    return _probe;
}

////////////////////////////////////////////////////////////////////

int ProbeFormBridge::numValues() const
{
    return _numValues;
}

////////////////////////////////////////////////////////////////////

string ProbeFormBridge::quantity() const
{
    return _quantity;
}

////////////////////////////////////////////////////////////////////

string ProbeFormBridge::projectedQuantity() const
{
    return _projectedQuantity;
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

void ProbeFormBridge::addColumnHeaders(TextOutFile* outfile) const
{
    if (!_addColumnHeaders) throw FATALERROR("Callback function addColumnHeaders has not been set");
    _addColumnHeaders(outfile);
}

////////////////////////////////////////////////////////////////////

string ProbeFormBridge::fileid() const
{
    return _fileid;
}

////////////////////////////////////////////////////////////////////

void ProbeFormBridge::valuesInCell(int m, Array& values) const
{
    if (!_valuesInCell) throw FATALERROR("Callback function valuesInCell has not been set");
    _valuesInCell(m, values);
}

////////////////////////////////////////////////////////////////////

void ProbeFormBridge::valuesAtPosition(Position bfr, Array& values) const
{
    if (_valuesAtPosition)
    {
        _valuesAtPosition(bfr, values);
    }
    else if (_valuesInCell)
    {
        int m = _grid ? _grid->cellIndex(bfr) : -1;
        if (m >= 0)
            _valuesInCell(m, values);
        else
            values = 0.;
    }
    else
    {
        throw FATALERROR("Callback functions valuesAtPosition nor valuesInCell have been set");
    }
}

////////////////////////////////////////////////////////////////////

void ProbeFormBridge::valuesAlongPath(Position bfr, Direction bfk, Array& values) const
{
    if (_valuesAlongPath)
    {
        _valuesAlongPath(bfr, bfk, values);
    }
    else if (_valuesInCell)
    {
        values = 0.;
        if (_grid)
        {
            // get a segment generator and initialize the path
            auto generator = _grid->createPathSegmentGenerator();
            SpatialGridPath path(bfr, bfk);

            // allocate temporary room for values of a single cell
            Array cellValues(_numValues);

            // accumulate values along the path
            if (_quantity != _projectedQuantity)
            {
                generator->start(&path);
                while (generator->next())
                {
                    int m = generator->m();
                    if (m >= 0)
                    {
                        _valuesInCell(m, cellValues);
                        values += cellValues * generator->ds();
                    }
                }
            }

            // average values along the path using weights provided by probe
            else
            {
                if (!_weightInCell) throw FATALERROR("Callback function weightInCell has not been set");

                double totalWeight = 0.;
                generator->start(&path);
                while (generator->next())
                {
                    int m = generator->m();
                    if (m >= 0)
                    {
                        _valuesInCell(m, cellValues);
                        double weight = _weightInCell(m) * generator->ds();
                        values += cellValues * weight;
                        totalWeight += weight;
                    }
                }
                values /= totalWeight;
            }
        }
    }
    else
    {
        throw FATALERROR("Callback functions valuesAlongPath nor valuesInCell have been set");
    }
}

////////////////////////////////////////////////////////////////////
