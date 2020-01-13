/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AbstractItemPropertyHandler.hpp"
#include "Item.hpp"
#include "NameManager.hpp"
#include "PropertyDef.hpp"
#include "SchemaDef.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

bool AbstractItemPropertyHandler::isValidValue(string value) const
{
    return !value.empty() && schema()->inherits(value, baseType())
           && nameManager()->evaluateBoolean(schema()->allowed(value));
}

////////////////////////////////////////////////////////////////////

bool AbstractItemPropertyHandler::isCompound() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

string AbstractItemPropertyHandler::baseType() const
{
    return property()->base();
}

////////////////////////////////////////////////////////////////////

string AbstractItemPropertyHandler::defaultType() const
{
    return nameManager()->evaluateConditionalValue(property()->defaultValue());
}

////////////////////////////////////////////////////////////////////

vector<string> AbstractItemPropertyHandler::allowedAndDisplayedDescendants()
{
    vector<string> descendants;
    for (auto candidate : schema()->descendants(property()->base()))
    {
        if (nameManager()->evaluateBoolean(schema()->allowedAndDisplayed(candidate))) descendants.push_back(candidate);
    }
    return descendants;
}

////////////////////////////////////////////////////////////////////
