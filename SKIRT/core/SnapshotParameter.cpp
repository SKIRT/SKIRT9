/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SnapshotParameter.hpp"

//////////////////////////////////////////////////////////////////////

SnapshotParameter::SnapshotParameter(Identifier identifier, string description, string quantity, string defaultUnit)
    : _identifier(identifier), _description(description), _quantity(quantity), _defaultUnit(defaultUnit)
{}

//////////////////////////////////////////////////////////////////////

SnapshotParameter SnapshotParameter::initialMass()
{
    return SnapshotParameter(Identifier::InitialMass, "initial mass", "mass", "Msun");
}

//////////////////////////////////////////////////////////////////////

SnapshotParameter SnapshotParameter::currentMass()
{
    return SnapshotParameter(Identifier::CurrentMass, "current mass", "mass", "Msun");
}

//////////////////////////////////////////////////////////////////////

SnapshotParameter SnapshotParameter::metallicity()
{
    return SnapshotParameter(Identifier::Metallicity, "metallicity", string(), string());
}

//////////////////////////////////////////////////////////////////////

SnapshotParameter SnapshotParameter::age()
{
    return SnapshotParameter(Identifier::Age, "age", "time", "yr");
}

//////////////////////////////////////////////////////////////////////

SnapshotParameter SnapshotParameter::temperature()
{
    return SnapshotParameter(Identifier::Temperature, "temperature", "temperature", "K");
}

//////////////////////////////////////////////////////////////////////

SnapshotParameter SnapshotParameter::custom(string description, string quantity, string defaultUnit)
{
    return SnapshotParameter(Identifier::Custom, description, quantity, defaultUnit);
}

//////////////////////////////////////////////////////////////////////
