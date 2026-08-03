/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Configuration.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

Configuration::Configuration(SimulationItem* parent)
{
    parent->addChild(this);
}

////////////////////////////////////////////////////////////////////

void Configuration::setEmulationMode()
{
    _emulationMode = true;
    _hasPrimaryIterations = false;
    _hasSecondaryIterations = false;
    _hasMergedIterations = false;
    _minPrimaryIterations = 0;
    _maxPrimaryIterations = 0;
    _minSecondaryIterations = 0;
    _maxSecondaryIterations = 0;
    _numPrimaryPackets = 0;
    _numPrimaryIterationPackets = 0;
    _numSecondaryPackets = 0;
    _numSecondaryIterationPackets = 0;
    _hasDynamicStateRecipes = false;
    _hasPrimaryDynamicStateMedia = false;
    _hasSecondaryDynamicStateMedia = false;
    _hasPrimaryDynamicState = false;
    _hasSecondaryDynamicState = false;
}

////////////////////////////////////////////////////////////////////

WavelengthGrid* Configuration::wavelengthGrid(WavelengthGrid* localWavelengthGrid) const
{
    auto result = localWavelengthGrid && !_oligochromatic ? localWavelengthGrid : _defaultWavelengthGrid;
    if (!result) throw FATALERROR("Cannot find a wavelength grid for instrument or probe");
    return result;
}

////////////////////////////////////////////////////////////////////
