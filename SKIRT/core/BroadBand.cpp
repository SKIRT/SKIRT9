/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BroadBand.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

BroadBand::BroadBand(SimulationItem* parent, string bandName)
{
    parent->addChild(this);
    _bandName = bandName;
    setup();
}

////////////////////////////////////////////////////////////////////

void BroadBand::setupSelfBefore()
{
    Band::setupSelfBefore();

    // locate the table (allow underscores, extra spaces, and lower/upper case)
    string bandName = StringUtils::toUpper(StringUtils::squeeze(StringUtils::replace(_bandName, "_", " ")));
    auto segments = StringUtils::split(bandName, " ");
    for (auto& segment : segments) segment += "_";
    string resourceName = FilePaths::resourceName("_BroadBand.stab", segments);

    // open the table
    _table.open(this, resourceName, "lambda(m)", "T(1)", false);
}

////////////////////////////////////////////////////////////////////

size_t BroadBand::dataSize() const
{
    return _table.size();
}

////////////////////////////////////////////////////////////////////

const double* BroadBand::wavelengthData() const
{
    return _table.axisData();
}

////////////////////////////////////////////////////////////////////

const double* BroadBand::transmissionData() const
{
    return _table.quantityData();
}

////////////////////////////////////////////////////////////////////
