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
#include "BlackBodySEDFamily.hpp"
#include "BoxClipGeometryDecorator.hpp"
#include "BrokenExpDiskGeometry.hpp"
#include "BruzualCharlotSED.hpp"
#include "BruzualCharlotSEDFamily.hpp"
#include "ClumpyGeometryDecorator.hpp"
#include "CombineGeometryDecorator.hpp"
#include "ConicalShellGeometry.hpp"
#include "CubicalBackgroundSource.hpp"
#include "CubicSplineSmoothingKernel.hpp"
#include "CylindricalClipGeometryDecorator.hpp"
#include "DefaultMediaDensityCutsProbe.hpp"
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
#include "ListWavelengthDistribution.hpp"
#include "ListWavelengthGrid.hpp"
#include "ListSED.hpp"
#include "LogWavelengthDistribution.hpp"
#include "LogWavelengthGrid.hpp"
#include "LuminosityProbe.hpp"
#include "MappingsSED.hpp"
#include "MappingsSEDFamily.hpp"
#include "MarastonSED.hpp"
#include "MarastonSEDFamily.hpp"
#include "MollweideProjection.hpp"
#include "MonteCarloSimulation.hpp"
#include "NestedLogWavelengthGrid.hpp"
#include "NetzerAngularDistribution.hpp"
#include "NoPolarizationProfile.hpp"
#include "OffsetGeometryDecorator.hpp"
#include "ParaboloidGeometry.hpp"
#include "ParaboloidShellGeometry.hpp"
#include "ParticleGeometry.hpp"
#include "ParticleSource.hpp"
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
#include "ScaledGaussianSmoothingKernel.hpp"
#include "SersicGeometry.hpp"
#include "ShellGeometry.hpp"
#include "SineSquarePolarizationProfile.hpp"
#include "SourceSystem.hpp"
#include "SpecificLuminosityNormalization.hpp"
#include "SphericalBackgroundSource.hpp"
#include "SphericalClipGeometryDecorator.hpp"
#include "SpheroidalGeometryDecorator.hpp"
#include "SpiralStructureGeometryDecorator.hpp"
#include "Starburst99SED.hpp"
#include "Starburst99SEDFamily.hpp"
#include "StellarSurfaceSource.hpp"
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
    ItemRegistry::add<ImportedSource>();
    ItemRegistry::add<ParticleSource>();
    ItemRegistry::add<CenteredSource>();
    ItemRegistry::add<StellarSurfaceSource>();
    ItemRegistry::add<CubicalBackgroundSource>();
    ItemRegistry::add<SphericalBackgroundSource>();

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
    ItemRegistry::add<FamilySED>();
    ItemRegistry::add<BruzualCharlotSED>();
    ItemRegistry::add<MarastonSED>();
    ItemRegistry::add<Starburst99SED>();
    ItemRegistry::add<MappingsSED>();
    ItemRegistry::add<TabulatedSED>();
    ItemRegistry::add<FileSED>();
    ItemRegistry::add<ListSED>();

    // SED families
    ItemRegistry::add<SEDFamily>();
    ItemRegistry::add<BlackBodySEDFamily>();
    ItemRegistry::add<BruzualCharlotSEDFamily>();
    ItemRegistry::add<MarastonSEDFamily>();
    ItemRegistry::add<Starburst99SEDFamily>();
    ItemRegistry::add<MappingsSEDFamily>();

    // Wavelength distributions
    ItemRegistry::add<WavelengthDistribution>();
    ItemRegistry::add<RangeWavelengthDistribution>();
    ItemRegistry::add<LinWavelengthDistribution>();
    ItemRegistry::add<LogWavelengthDistribution>();
    ItemRegistry::add<TabulatedWavelengthDistribution>();
    ItemRegistry::add<FileWavelengthDistribution>();
    ItemRegistry::add<ListWavelengthDistribution>();

    // Angular distributions
    ItemRegistry::add<AngularDistribution>();
    ItemRegistry::add<IsotropicAngularDistribution>();
    ItemRegistry::add<AxAngularDistribution>();
    ItemRegistry::add<LaserAngularDistribution>();
    ItemRegistry::add<NetzerAngularDistribution>();

    // Polarization profiles
    ItemRegistry::add<PolarizationProfile>();
    ItemRegistry::add<NoPolarizationProfile>();
    ItemRegistry::add<SineSquarePolarizationProfile>();

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
    ItemRegistry::add<UniformBoxGeometry>();
    ItemRegistry::add<ReadFitsGeometry>();
    ItemRegistry::add<ImportedGeometry>();
    ItemRegistry::add<ParticleGeometry>();

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
    ItemRegistry::add<CubicSplineSmoothingKernel>();
    ItemRegistry::add<ScaledGaussianSmoothingKernel>();
    ItemRegistry::add<UniformSmoothingKernel>();

    // wavelength grids
    ItemRegistry::add<WavelengthGrid>();
    ItemRegistry::add<LogWavelengthGrid>();
    ItemRegistry::add<NestedLogWavelengthGrid>();
    ItemRegistry::add<FileWavelengthGrid>();
    ItemRegistry::add<ListWavelengthGrid>();

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

    // probe system and probes
    ItemRegistry::add<ProbeSystem>();
    ItemRegistry::add<Probe>();
    ItemRegistry::add<WavelengthGridProbe>();
    ItemRegistry::add<LuminosityProbe>();
    ItemRegistry::add<LaunchedPacketsProbe>();
    ItemRegistry::add<DefaultMediaDensityCutsProbe>();

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

