/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SimulationItemRegistry.hpp"
#include "ItemRegistry.hpp"
#include "SkirtUnitDef.hpp"

// ---> add new items below in alphabetical order

#include "AllSkyInstrument.hpp"
#include "ExtragalacticUnits.hpp"
#include "FileWavelengthGrid.hpp"
#include "FrameInstrument.hpp"
#include "FullInstrument.hpp"
#include "InstrumentSystem.hpp"
#include "ListWavelengthGrid.hpp"
#include "LogWavelengthGrid.hpp"
#include "MonteCarloSimulation.hpp"
#include "NestedLogWavelengthGrid.hpp"
#include "PerspectiveInstrument.hpp"
#include "Random.hpp"
#include "SEDInstrument.hpp"
#include "SIUnits.hpp"
#include "StellarUnits.hpp"

////////////////////////////////////////////////////////////////////

SimulationItemRegistry::SimulationItemRegistry(string version, string format)
{
    // start a new schema
    ItemRegistry::beginSchema("SKIRT", "a SKIRT parameter file", version, "ski",
                              "skirt-simulation-hierarchy", "MonteCarloSimulation", format,
                              "http://www.skirt.ugent.be/skirt");

    // add the SKIRT unit definitions
    ItemRegistry::addUnitDef<SkirtUnitDef>();

    // add the SKIRT simulation items
    ItemRegistry::add<SimulationItem>();

    // ---> add new items in the order you want them to appear in choice lists for the user

    // basic building blocks
    ItemRegistry::add<Simulation>();
    ItemRegistry::add<Random>();
    ItemRegistry::add<Units>();
    ItemRegistry::add<SIUnits>();
    ItemRegistry::add<StellarUnits>();
    ItemRegistry::add<ExtragalacticUnits>();

    // wavelength grids
    ItemRegistry::add<WavelengthGrid>();
    ItemRegistry::add<ListWavelengthGrid>();
    ItemRegistry::add<LogWavelengthGrid>();
    ItemRegistry::add<NestedLogWavelengthGrid>();
    ItemRegistry::add<FileWavelengthGrid>();

    // instrument system and instruments
    ItemRegistry::add<InstrumentSystem>();
    ItemRegistry::add<Instrument>();
    ItemRegistry::add<DistantInstrument>();
    ItemRegistry::add<SEDInstrument>();
    ItemRegistry::add<FrameInstrument>();
    ItemRegistry::add<FullInstrument>();
    ItemRegistry::add<AllSkyInstrument>();
    ItemRegistry::add<PerspectiveInstrument>();

    // Monte Carlo simulations
    ItemRegistry::add<MonteCarloSimulation>();
}

////////////////////////////////////////////////////////////////////

const SchemaDef* SimulationItemRegistry::getSchemaDef()
{
    return ItemRegistry::getSchemaDef("SKIRT");
}

////////////////////////////////////////////////////////////////////

SimulationItemRegistry::~SimulationItemRegistry()
{
    ItemRegistry::finalize();
}

////////////////////////////////////////////////////////////////////

