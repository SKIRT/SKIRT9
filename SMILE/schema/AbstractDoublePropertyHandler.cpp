/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AbstractDoublePropertyHandler.hpp"
#include "EnumPropertyHandler.hpp"
#include "FatalError.hpp"
#include "Item.hpp"
#include "PropertyDef.hpp"
#include "SchemaDef.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

double AbstractDoublePropertyHandler::minValue() const
{
    // get minimum value as a non-empty string
    string value = StringUtils::squeeze(property()->minValue());
    if (value.empty()) return -std::numeric_limits<double>::infinity();

    // remove leading "[" or "]"
    if (value.front() == '[' || value.front() == ']') value.erase(0, 1);

    // convert to double
    if (isValidDouble(value)) return toDouble(value);
    throw FATALERROR("Invalid minimum value '" + value + "' for property '" + property()->name() + "' of type '"
                     + target()->type() + "'");
}

////////////////////////////////////////////////////////////////////

double AbstractDoublePropertyHandler::maxValue() const
{
    // get maximum value as a non-empty string
    string value = StringUtils::squeeze(property()->maxValue());
    if (value.empty()) return std::numeric_limits<double>::infinity();

    // remove trailing "[" or "]"
    if (value.back() == '[' || value.back() == ']') value.pop_back();

    // convert to double
    if (isValidDouble(value)) return toDouble(value);
    throw FATALERROR("Invalid maximum value '" + value + "' for property '" + property()->name() + "' of type '"
                     + target()->type() + "'");
}

////////////////////////////////////////////////////////////////////

bool AbstractDoublePropertyHandler::minOpen() const
{
    string value = StringUtils::squeeze(property()->minValue());
    return value.empty() || value.front() == ']';
}

////////////////////////////////////////////////////////////////////

bool AbstractDoublePropertyHandler::maxOpen() const
{
    string value = StringUtils::squeeze(property()->maxValue());
    return value.empty() || value.back() == '[';
}

////////////////////////////////////////////////////////////////////

string AbstractDoublePropertyHandler::rangeDescription() const
{
    return (minOpen() ? "]" : "[") + toString(minValue()) + "," + toString(maxValue()) + (maxOpen() ? "[" : "]");
}

////////////////////////////////////////////////////////////////////

bool AbstractDoublePropertyHandler::isInRange(double value) const
{
    if (!std::isfinite(value)) return false;

    if (minOpen())
    {
        if (value <= minValue()) return false;
    }
    else
    {
        if (value < minValue()) return false;
    }

    if (maxOpen())
    {
        if (value >= maxValue()) return false;
    }
    else
    {
        if (value > maxValue()) return false;
    }

    return true;
}

////////////////////////////////////////////////////////////////////

bool AbstractDoublePropertyHandler::isInRange(const vector<double>& value) const
{
    return std::all_of(value.cbegin(), value.cend(), [this](double x) { return isInRange(x); });
}

////////////////////////////////////////////////////////////////////

bool AbstractDoublePropertyHandler::isValidDouble(string value) const
{
    // split into segments; must have exactly one or two segments
    vector<string> segments = StringUtils::split(StringUtils::squeeze(value), " ");
    if (segments.empty() || segments.size() > 2) return false;

    // verify that the first segment is a valid number
    if (!StringUtils::isValidDouble(segments[0])) return false;

    // handle the physical quantity
    string qty = quantity();
    if (qty.empty())
    {
        // a dimensionless quantity may have no unit spec
        return segments.size() == 1;
    }
    else
    {
        // use the unit spec in the input value if there is one; otherwise use the unit system to get a default unit
        string unit = segments.size() == 2 ? segments[1] : unitSystem();

        // verify that the combination works
        // since this is probably an input error rather than a programming error, we don't throw an exception
        return schema()->has(qty, unit);
    }
}

////////////////////////////////////////////////////////////////////

double AbstractDoublePropertyHandler::toDouble(string value) const
{
    // ensure that the string is a valid representation
    if (!isValidDouble(value)) return 0.;

    // split into one or two segments
    vector<string> segments = StringUtils::split(StringUtils::squeeze(value), " ");

    // convert the first segment to a number
    double result = StringUtils::toDouble(segments[0]);

    // handle the physical quantity, if needed
    string qty = quantity();
    if (!qty.empty())
    {
        // use the unit spec in the input value if there is one; otherwise use the unit system to get a default unit
        string unit = segments.size() == 2 ? segments[1] : unitSystem();
        result = schema()->in(qty, unit, result);
    }

    return result;
}

////////////////////////////////////////////////////////////////////

string AbstractDoublePropertyHandler::toString(double value) const
{
    // get unit specification and convert value to external units if needed
    string unit;
    string qty = quantity();
    if (!qty.empty())
    {
        string system = unitSystem();
        value = schema()->out(qty, system, value);  // overwrite incoming value
        unit = " " + schema()->unit(qty, system);   // include separating space
    }

    return StringUtils::toString(value) + unit;
}

////////////////////////////////////////////////////////////////////

bool AbstractDoublePropertyHandler::isValidDoubleList(string value) const
{
    // split into comma-seperated segments and ensure there is at least one segment
    vector<string> segments = StringUtils::split(value, ",");
    if (segments.empty()) return false;

    // validate each segment
    for (string segment : segments)
        if (!isValidDouble(segment)) return false;
    return true;
}

////////////////////////////////////////////////////////////////////

vector<double> AbstractDoublePropertyHandler::toDoubleList(string value) const
{
    vector<double> result;
    if (isValidDoubleList(value))
    {
        for (string segment : StringUtils::split(value, ",")) result.push_back(toDouble(segment));
    }
    return result;
}

////////////////////////////////////////////////////////////////////

string AbstractDoublePropertyHandler::toString(vector<double> value) const
{
    string result;
    string separator;
    for (double number : value)
    {
        result += separator + toString(number);
        separator = ", ";
    }
    return result;
}

////////////////////////////////////////////////////////////////////

string AbstractDoublePropertyHandler::quantity() const
{
    string result = property()->quantity();

    // if the attribute value starts with an at sign,
    // the quantity is determined as the value of the indicated enumeration property
    if (StringUtils::startsWith(result, "@"))
    {
        // construct a handler for the target property and get its enumeration value
        auto handler = schema()->createPropertyHandler(target(), result.substr(1), nameManager());
        EnumPropertyHandler* enumHandler = dynamic_cast<EnumPropertyHandler*>(handler.get());

        // if the property has the wrong type, or its value is an unknown quantity,
        // return the empty string, which means "dimensionless"
        if (enumHandler && schema()->has(enumHandler->value()))
            result = enumHandler->value();
        else
            result.clear();
    }
    return result;
}

////////////////////////////////////////////////////////////////////

string AbstractDoublePropertyHandler::unitSystem() const
{
    // get the base type; or if there is only one unit system, its name
    string baseType = schema()->unitSystemBase();
    if (!schema()->hasMultipleUnitSystems()) return baseType;

    // loop over all ancestors
    Item* ancestor = target();
    while (ancestor)
    {
        // test the ancestor
        if (schema()->inherits(ancestor->type(), baseType)) return ancestor->type();

        // test its children
        for (auto child : ancestor->children())
        {
            if (schema()->inherits(child->type(), baseType)) return child->type();
        }

        // next ancestor
        ancestor = ancestor->parent();
    }
    throw FATALERROR("No item of type " + baseType + " found in hierarchy");
}

////////////////////////////////////////////////////////////////////
