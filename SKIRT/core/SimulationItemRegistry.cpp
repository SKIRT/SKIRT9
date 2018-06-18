/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SimulationItemRegistry.hpp"
#include "ItemRegistry.hpp"
#include "SkirtUnitDef.hpp"

// ---> add new items below in alphabetical order

#include "AllSkyInstrument.hpp"
#include "BlackBodySED.hpp"
#include "BoxClipGeometryDecorator.hpp"
#include "BrokenExpDiskGeometry.hpp"
#include "ClumpyGeometryDecorator.hpp"
#include "CombineGeometryDecorator.hpp"
#include "ConicalShellGeometry.hpp"
#include "CubicSplineSmoothingKernel.hpp"
#include "CylindricalClipGeometryDecorator.hpp"
#include "EinastoGeometry.hpp"
#include "ExpDiskGeometry.hpp"
#include "ExtragalacticUnits.hpp"
#include "FileSED.hpp"
#include "FileWavelengthDistribution.hpp"
#include "FileWavelengthGrid.hpp"
#include "FrameInstrument.hpp"
#include "FullInstrument.hpp"
#include "GammaGeometry.hpp"
#include "GaussianGeometry.hpp"
#include "GeometricSource.hpp"
#include "HammerAitoffProjection.hpp"
#include "HyperboloidGeometry.hpp"
#include "HyperboloidShellGeometry.hpp"
#include "InstrumentSystem.hpp"
#include "IntegratedLuminosityNormalization.hpp"
#include "IsotropicAngularDistribution.hpp"
#include "LaserAngularDistribution.hpp"
#include "LaunchedPacketsProbe.hpp"
#include "LinWavelengthDistribution.hpp"
#include "ListWavelengthGrid.hpp"
#include "LogWavelengthDistribution.hpp"
#include "LogWavelengthGrid.hpp"
#include "LuminosityProbe.hpp"
#include "MollweideProjection.hpp"
#include "MonteCarloSimulation.hpp"
#include "NestedLogWavelengthGrid.hpp"
#include "NetzerAngularDistribution.hpp"
#include "OffsetGeometryDecorator.hpp"
#include "ParaboloidGeometry.hpp"
#include "ParaboloidShellGeometry.hpp"
#include "PerspectiveInstrument.hpp"
#include "PlummerGeometry.hpp"
#include "PointSource.hpp"
#include "ProbeSystem.hpp"
#include "PseudoSersicGeometry.hpp"
#include "QuasarSED.hpp"
#include "Random.hpp"
#include "ReadFitsGeometry.hpp"
#include "RingGeometry.hpp"
#include "RotateGeometryDecorator.hpp"
#include "SEDInstrument.hpp"
#include "SIUnits.hpp"
#include "SersicGeometry.hpp"
#include "ShellGeometry.hpp"
#include "SourceSystem.hpp"
#include "SpecificLuminosityNormalization.hpp"
#include "SphericalClipGeometryDecorator.hpp"
#include "SpheroidalGeometryDecorator.hpp"
#include "SpiralStructureGeometryDecorator.hpp"
#include "StellarUnits.hpp"
#include "SunSED.hpp"
#include "TorusGeometry.hpp"
#include "TriaxialGeometryDecorator.hpp"
#include "UniformBoxGeometry.hpp"
#include "UniformSmoothingKernel.hpp"
#include "WavelengthGridProbe.hpp"

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

    // source system and sources
    ItemRegistry::add<SourceSystem>();
    ItemRegistry::add<Source>();
    ItemRegistry::add<NormalizedSource>();
    ItemRegistry::add<PointSource>();
    ItemRegistry::add<GeometricSource>();

    // luminosity normalizations
    ItemRegistry::add<LuminosityNormalization>();
    ItemRegistry::add<IntegratedLuminosityNormalization>();
    ItemRegistry::add<SpecificLuminosityNormalization>();

    // SEDs
    ItemRegistry::add<SED>();
    ItemRegistry::add<BlackBodySED>();
    ItemRegistry::add<ResourceSED>();
    ItemRegistry::add<SunSED>();
    ItemRegistry::add<QuasarSED>();
    ItemRegistry::add<FileSED>();

    // Wavelength distributions
    ItemRegistry::add<WavelengthDistribution>();
    ItemRegistry::add<RangeWavelengthDistribution>();
    ItemRegistry::add<LinWavelengthDistribution>();
    ItemRegistry::add<LogWavelengthDistribution>();
    ItemRegistry::add<FileWavelengthDistribution>();

    // Angular distributions
    ItemRegistry::add<AngularDistribution>();
    ItemRegistry::add<IsotropicAngularDistribution>();
    ItemRegistry::add<LaserAngularDistribution>();
    ItemRegistry::add<NetzerAngularDistribution>();

    // geometries
    ItemRegistry::add<Geometry>();
    ItemRegistry::add<SpheGeometry>();
    ItemRegistry::add<PlummerGeometry>();
    ItemRegistry::add<GammaGeometry>();
    ItemRegistry::add<SersicGeometry>();
    ItemRegistry::add<PseudoSersicGeometry>();
    ItemRegistry::add<EinastoGeometry>();
    ItemRegistry::add<GaussianGeometry>();
    ItemRegistry::add<ShellGeometry>();
    ItemRegistry::add<AxGeometry>();
    ItemRegistry::add<SepAxGeometry>();
    ItemRegistry::add<ExpDiskGeometry>();
    ItemRegistry::add<BrokenExpDiskGeometry>();
    ItemRegistry::add<RingGeometry>();
    ItemRegistry::add<TorusGeometry>();
    ItemRegistry::add<ConicalShellGeometry>();
    ItemRegistry::add<ParaboloidGeometry>();
    ItemRegistry::add<ParaboloidShellGeometry>();
    ItemRegistry::add<HyperboloidGeometry>();
    ItemRegistry::add<HyperboloidShellGeometry>();
    ItemRegistry::add<GenGeometry>();
    ItemRegistry::add<BoxGeometry>();
    ItemRegistry::add<UniformBoxGeometry>();
    ItemRegistry::add<ReadFitsGeometry>();

    // geometry decorators
    ItemRegistry::add<OffsetGeometryDecorator>();
    ItemRegistry::add<RotateGeometryDecorator>();
    ItemRegistry::add<SpheroidalGeometryDecorator>();
    ItemRegistry::add<TriaxialGeometryDecorator>();
    ItemRegistry::add<ClipGeometryDecorator>();
    ItemRegistry::add<SphericalClipGeometryDecorator>();
    ItemRegistry::add<CylindricalClipGeometryDecorator>();
    ItemRegistry::add<BoxClipGeometryDecorator>();
    ItemRegistry::add<SpiralStructureGeometryDecorator>();
    ItemRegistry::add<ClumpyGeometryDecorator>();
    ItemRegistry::add<CombineGeometryDecorator>();

    // smoothing kernels
    ItemRegistry::add<SmoothingKernel>();
    ItemRegistry::add<UniformSmoothingKernel>();
    ItemRegistry::add<CubicSplineSmoothingKernel>();

    // wavelength grids
    ItemRegistry::add<WavelengthGrid>();
    ItemRegistry::add<ListWavelengthGrid>();
    ItemRegistry::add<LogWavelengthGrid>();
    ItemRegistry::add<NestedLogWavelengthGrid>();
    ItemRegistry::add<FileWavelengthGrid>();

    // probe system and probes
    ItemRegistry::add<ProbeSystem>();
    ItemRegistry::add<Probe>();
    ItemRegistry::add<WavelengthGridProbe>();
    ItemRegistry::add<LuminosityProbe>();
    ItemRegistry::add<LaunchedPacketsProbe>();

    // instrument system and instruments
    ItemRegistry::add<InstrumentSystem>();
    ItemRegistry::add<Instrument>();
    ItemRegistry::add<DistantInstrument>();
    ItemRegistry::add<SEDInstrument>();
    ItemRegistry::add<FrameInstrument>();
    ItemRegistry::add<FullInstrument>();
    ItemRegistry::add<AllSkyInstrument>();
    ItemRegistry::add<PerspectiveInstrument>();

    // all-sky projections
    ItemRegistry::add<AllSkyProjection>();
    ItemRegistry::add<HammerAitoffProjection>();
    ItemRegistry::add<MollweideProjection>();

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

