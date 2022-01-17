/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "UnitDef.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

UnitDef::UnitDef() {}

////////////////////////////////////////////////////////////////////

void UnitDef::addUnit(string quantity, string unit, double factor, double power, double offset)
{
    _quantities[quantity][unit] = std::make_tuple(factor, power, offset);
}

////////////////////////////////////////////////////////////////////

void UnitDef::addDefaultUnit(string unitSystem, string quantity, string unit)
{
    _unitSystems[unitSystem][quantity] = unit;
}

////////////////////////////////////////////////////////////////////

bool UnitDef::has(string qty) const
{
    return _quantities.count(qty) != 0;
}

////////////////////////////////////////////////////////////////////

bool UnitDef::has(string qty, string unit) const
{
    // if the unit argument represents a unit system, replace it by the default unit for the quantity
    if (_unitSystems.count(unit) && _unitSystems.at(unit).count(qty)) unit = _unitSystems.at(unit).at(qty);

    // check whether the unit is defined for the quantity
    return _quantities.count(qty) && _quantities.at(qty).count(unit);
}

////////////////////////////////////////////////////////////////////

std::tuple<double, double, double> UnitDef::def(string qty, string unit) const
{
    // if the unit argument represents a unit system, replace it by the default unit for the quantity
    if (_unitSystems.count(unit) && _unitSystems.at(unit).count(qty)) unit = _unitSystems.at(unit).at(qty);

    // if the unit is defined for the quantity, return its definition
    if (_quantities.count(qty) && _quantities.at(qty).count(unit))
    {
        return _quantities.at(qty).at(unit);
    }

    // otherwise report the error
    throw FATALERROR("Unknow quantity " + qty + " and/or unit (system) " + unit);
}

////////////////////////////////////////////////////////////////////

double UnitDef::in(string qty, string unit, double value) const
{
    // if the unit argument represents a unit system, replace it by the default unit for the quantity
    if (_unitSystems.count(unit) && _unitSystems.at(unit).count(qty)) unit = _unitSystems.at(unit).at(qty);

    // if the unit is defined for the quantity, perform the conversion
    if (_quantities.count(qty) && _quantities.at(qty).count(unit))
    {
        double factor, power, offset;
        std::tie(factor, power, offset) = _quantities.at(qty).at(unit);
        if (power != 1.) value = pow(value, power);
        return factor * value + offset;
    }

    // otherwise report the error
    throw FATALERROR("Unknow quantity " + qty + " and/or unit (system) " + unit);
}

////////////////////////////////////////////////////////////////////

double UnitDef::out(string qty, string unit, double value) const
{
    // if the unit argument represents a unit system, replace it by the default unit for the quantity
    if (_unitSystems.count(unit) && _unitSystems.at(unit).count(qty)) unit = _unitSystems.at(unit).at(qty);

    // if the unit is defined for the quantity, perform the conversion
    if (_quantities.count(qty) && _quantities.at(qty).count(unit))
    {
        double factor, power, offset;
        std::tie(factor, power, offset) = _quantities.at(qty).at(unit);
        value = (value - offset) / factor;
        if (power != 1.) value = pow(value, 1. / power);
        return value;
    }

    // otherwise report the error
    throw FATALERROR("Unknow quantity " + qty + " and/or unit (system) " + unit);
}

////////////////////////////////////////////////////////////////////

string UnitDef::unit(string qty, string unitSystem) const
{
    if (_unitSystems.count(unitSystem) && _unitSystems.at(unitSystem).count(qty))
        return _unitSystems.at(unitSystem).at(qty);

    throw FATALERROR("Unknow quantity " + qty + " and/or unit system " + unitSystem);
}

////////////////////////////////////////////////////////////////////
