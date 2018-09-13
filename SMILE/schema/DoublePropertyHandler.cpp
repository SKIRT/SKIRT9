/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DoublePropertyHandler.hpp"
#include "Item.hpp"
#include "NameManager.hpp"
#include "PropertyDef.hpp"
#include "PropertyHandlerVisitor.hpp"

////////////////////////////////////////////////////////////////////

bool DoublePropertyHandler::isValidValue(string value) const
{
    return isValidDouble(value);
}

////////////////////////////////////////////////////////////////////

void DoublePropertyHandler::insertNames()
{
    if (value() != 0.) nameManager()->insert(property()->name());
    nameManager()->insertFromConditionalValue(property()->insert());
}

////////////////////////////////////////////////////////////////////

void DoublePropertyHandler::acceptVisitor(PropertyHandlerVisitor* visitor)
{
    visitor->visitPropertyHandler(this);
}

////////////////////////////////////////////////////////////////////

double DoublePropertyHandler::defaultValue() const
{
    return toDouble(nameManager()->evaluateConditionalValue(property()->defaultValue()));
}

////////////////////////////////////////////////////////////////////

double DoublePropertyHandler::value() const
{
    return target()->getDoubleProperty(property());
}

////////////////////////////////////////////////////////////////////

void DoublePropertyHandler::setValue(double value)
{
    target()->setDoubleProperty(property(), value);
    setChanged();
}

////////////////////////////////////////////////////////////////////
