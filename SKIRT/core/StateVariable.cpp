/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StateVariable.hpp"

//////////////////////////////////////////////////////////////////////

StateVariable::StateVariable(Identifier identifier, int customIndex, string description, string quantity, char format)
    : _identifier(identifier), _customIndex(customIndex), _description(description), _quantity(quantity),
      _format(format)
{}

//////////////////////////////////////////////////////////////////////

StateVariable StateVariable::volume()
{
    return StateVariable(Identifier::Volume, 0, "volume", "volume", 'e');
}

//////////////////////////////////////////////////////////////////////

StateVariable StateVariable::bulkVelocity()
{
    return StateVariable(Identifier::BulkVelocity, 0, "bulk velocity", "velocity", 'e');
}

//////////////////////////////////////////////////////////////////////

StateVariable StateVariable::magneticField()
{
    return StateVariable(Identifier::MagneticField, 0, "magnetic field", "magneticfield", 'e');
}

//////////////////////////////////////////////////////////////////////

StateVariable StateVariable::numberDensity()
{
    return StateVariable(Identifier::NumberDensity, 0, "number density", "numbervolumedensity", 'e');
}

//////////////////////////////////////////////////////////////////////

StateVariable StateVariable::metallicity()
{
    return StateVariable(Identifier::Metallicity, 0, "metallicity", "", 'e');
}

//////////////////////////////////////////////////////////////////////

StateVariable StateVariable::temperature()
{
    return StateVariable(Identifier::Temperature, 0, "temperature", "temperature", 'e');
}

//////////////////////////////////////////////////////////////////////

StateVariable StateVariable::custom(int customIndex, string description, string quantity, char format)
{
    return StateVariable(Identifier::Custom, customIndex, description, quantity, format);
}

//////////////////////////////////////////////////////////////////////
