/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "IntPropertyHandler.hpp"
#include "Item.hpp"
#include "PropertyDef.hpp"
#include "PropertyHandlerVisitor.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // use a nice "round" maximum close to 2**31
    constexpr int MAXINT = 2*1000*1000*1000;
}

////////////////////////////////////////////////////////////////////

bool IntPropertyHandler::isTrueInCondition() const
{
    return value() != 0;
}

////////////////////////////////////////////////////////////////////

bool IntPropertyHandler::hasDefaultValue() const
{
    return StringUtils::isValidInt(property()->defaultValue());
}

////////////////////////////////////////////////////////////////////

void IntPropertyHandler::acceptVisitor(PropertyHandlerVisitor* visitor)
{
    visitor->visitPropertyHandler(this);
}

////////////////////////////////////////////////////////////////////

int IntPropertyHandler::defaultValue() const
{
    return StringUtils::toInt(property()->defaultValue());
}

////////////////////////////////////////////////////////////////////

int IntPropertyHandler::minValue() const
{
    string value = property()->minValue();
    if (StringUtils::isValidInt(value)) return StringUtils::toInt(value);
    else return -MAXINT;
}

////////////////////////////////////////////////////////////////////

int IntPropertyHandler::maxValue() const
{
    string value = property()->maxValue();
    if (StringUtils::isValidInt(value)) return StringUtils::toInt(value);
    else return MAXINT;
}

////////////////////////////////////////////////////////////////////

int IntPropertyHandler::value() const
{
    return target()->getIntProperty(property());
}

////////////////////////////////////////////////////////////////////

void IntPropertyHandler::setValue(int value)
{
    target()->setIntProperty(property(), value);
    setChanged();
}

////////////////////////////////////////////////////////////////////
