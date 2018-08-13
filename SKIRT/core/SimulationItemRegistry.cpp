/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SimulationItemRegistry.hpp"
#include "ItemRegistry.hpp"
#include "SkirtUnitDef.hpp"

// ---> add new items below in alphabetical order

#include "AdaptiveMeshGeometry.hpp"
#include "AdaptiveMeshSource.hpp"
#include "AllSkyInstrument.hpp"
#include "BandLuminosityNormalization.hpp"
#include "BandWavelengthGrid.hpp"
#include "BlackBodySED.hpp"
#include "BlackBodySEDFamily.hpp"
#include "BoxClipGeometryDecorator.hpp"
#include "BroadBand.hpp"
#include "BrokenExpDiskGeometry.hpp"
#include "BruzualCharlotSED.hpp"
#include "BruzualCharlotSEDFamily.hpp"
#include "CartesianSpatialGrid.hpp"
#include "CastelliKuruczSED.hpp"
#include "CastelliKuruczSEDFamily.hpp"
#include "ClumpyGeometryDecorator.hpp"
#include "CombineGeometryDecorator.hpp"
#include "ConicalAngularDistribution.hpp"
#include "ConicalShellGeometry.hpp"
#include "CubicSplineSmoothingKernel.hpp"
#include "CubicalBackgroundSource.hpp"
#include "Cylinder2DSpatialGrid.hpp"
#include "CylindricalClipGeometryDecorator.hpp"
#include "DefaultMediaDensityCutsProbe.hpp"
#include "EinastoGeometry.hpp"
#include "ElectronMix.hpp"
#include "ExpDiskGeometry.hpp"
#include "ExtinctionOnlyMode.hpp"
#include "ExtragalacticUnits.hpp"
#include "FileBand.hpp"
#include "FileMesh.hpp"
#include "FileSED.hpp"
#include "FileWavelengthDistribution.hpp"
#include "FileWavelengthGrid.hpp"
#include "FrameInstrument.hpp"
#include "FullInstrument.hpp"
#include "GammaGeometry.hpp"
#include "GaussianGeometry.hpp"
#include "GeometricMedium.hpp"
#include "GeometricSource.hpp"
#include "HammerAitoffProjection.hpp"
#include "HyperboloidGeometry.hpp"
#include "HyperboloidShellGeometry.hpp"
#include "InstrumentSystem.hpp"
#include "IntegratedLuminosityNormalization.hpp"
#include "IsotropicAngularDistribution.hpp"
#include "LaserAngularDistribution.hpp"
#include "LaunchedPacketsProbe.hpp"
#include "LinMesh.hpp"
#include "LinWavelengthDistribution.hpp"
#include "ListBand.hpp"
#include "ListSED.hpp"
#include "ListWavelengthDistribution.hpp"
#include "ListWavelengthGrid.hpp"
#include "LogMesh.hpp"
#include "LogWavelengthDistribution.hpp"
#include "LogWavelengthGrid.hpp"
#include "LuminosityProbe.hpp"
#include "MappingsSED.hpp"
#include "MappingsSEDFamily.hpp"
#include "MarastonSED.hpp"
#include "MarastonSEDFamily.hpp"
#include "MassColumnMaterialNormalization.hpp"
#include "MassMaterialNormalization.hpp"
#include "MeanTrustBenchmarkDustMix.hpp"
#include "MediumSystem.hpp"
#include "MollweideProjection.hpp"
#include "MonteCarloSimulation.hpp"
#include "MultiGaussianExpansionGeometry.hpp"
#include "NestedLogWavelengthGrid.hpp"
#include "NetzerAngularDistribution.hpp"
#include "NoPolarizationProfile.hpp"
#include "NumberColumnMaterialNormalization.hpp"
#include "OffsetGeometryDecorator.hpp"
#include "OpticalDepthMaterialNormalization.hpp"
#include "ParaboloidGeometry.hpp"
#include "ParaboloidShellGeometry.hpp"
#include "ParticleGeometry.hpp"
#include "ParticleSource.hpp"
#include "PerspectiveInstrument.hpp"
#include "PlummerGeometry.hpp"
#include "PointSource.hpp"
#include "PowMesh.hpp"
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
#include "SpatialGrid.hpp"
#include "SpecificLuminosityNormalization.hpp"
#include "Sphere1DSpatialGrid.hpp"
#include "Sphere2DSpatialGrid.hpp"
#include "SphericalBackgroundSource.hpp"
#include "SphericalClipGeometryDecorator.hpp"
#include "SpheroidalGeometryDecorator.hpp"
#include "SpiralStructureGeometryDecorator.hpp"
#include "Starburst99SED.hpp"
#include "Starburst99SEDFamily.hpp"
#include "StellarSurfaceSource.hpp"
#include "StellarUnits.hpp"
#include "SunSED.hpp"
#include "SymPowMesh.hpp"
#include "TTauriDiskGeometry.hpp"
#include "TorusGeometry.hpp"
#include "TriaxialGeometryDecorator.hpp"
#include "UniformBoxGeometry.hpp"
#include "UniformSmoothingKernel.hpp"
#include "VoronoiMeshGeometry.hpp"
#include "VoronoiMeshSource.hpp"
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

    // simulation modes
    ItemRegistry::add<SimulationMode>();
    ItemRegistry::add<NoMediaMode>();
    ItemRegistry::add<ExtinctionOnlyMode>();

    // source system and sources
    ItemRegistry::add<SourceSystem>();
    ItemRegistry::add<Source>();
    ItemRegistry::add<NormalizedSource>();
    ItemRegistry::add<PointSource>();
    ItemRegistry::add<GeometricSource>();
    ItemRegistry::add<ImportedSource>();
    ItemRegistry::add<ParticleSource>();
    ItemRegistry::add<MeshSource>();
    ItemRegistry::add<AdaptiveMeshSource>();
    ItemRegistry::add<VoronoiMeshSource>();
    ItemRegistry::add<CenteredSource>();
    ItemRegistry::add<StellarSurfaceSource>();
    ItemRegistry::add<CubicalBackgroundSource>();
    ItemRegistry::add<SphericalBackgroundSource>();

    // luminosity normalizations
    ItemRegistry::add<LuminosityNormalization>();
    ItemRegistry::add<IntegratedLuminosityNormalization>();
    ItemRegistry::add<SpecificLuminosityNormalization>();
    ItemRegistry::add<BandLuminosityNormalization>();

    // SEDs
    ItemRegistry::add<SED>();
    ItemRegistry::add<BlackBodySED>();
    ItemRegistry::add<ResourceSED>();
    ItemRegistry::add<SunSED>();
    ItemRegistry::add<QuasarSED>();
    ItemRegistry::add<FamilySED>();
    ItemRegistry::add<CastelliKuruczSED>();
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
    ItemRegistry::add<CastelliKuruczSEDFamily>();
    ItemRegistry::add<BruzualCharlotSEDFamily>();
    ItemRegistry::add<MarastonSEDFamily>();
    ItemRegistry::add<Starburst99SEDFamily>();
    ItemRegistry::add<MappingsSEDFamily>();

    // wavelength distributions
    ItemRegistry::add<WavelengthDistribution>();
    ItemRegistry::add<RangeWavelengthDistribution>();
    ItemRegistry::add<LinWavelengthDistribution>();
    ItemRegistry::add<LogWavelengthDistribution>();
    ItemRegistry::add<TabulatedWavelengthDistribution>();
    ItemRegistry::add<FileWavelengthDistribution>();
    ItemRegistry::add<ListWavelengthDistribution>();

    // bands
    ItemRegistry::add<Band>();
    ItemRegistry::add<BroadBand>();
    ItemRegistry::add<FileBand>();
    ItemRegistry::add<ListBand>();

    // angular distributions
    ItemRegistry::add<AngularDistribution>();
    ItemRegistry::add<IsotropicAngularDistribution>();
    ItemRegistry::add<AxAngularDistribution>();
    ItemRegistry::add<LaserAngularDistribution>();
    ItemRegistry::add<ConicalAngularDistribution>();
    ItemRegistry::add<NetzerAngularDistribution>();

    // polarization profiles
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
    ItemRegistry::add<TTauriDiskGeometry>();
    ItemRegistry::add<ConicalShellGeometry>();
    ItemRegistry::add<ParaboloidGeometry>();
    ItemRegistry::add<ParaboloidShellGeometry>();
    ItemRegistry::add<HyperboloidGeometry>();
    ItemRegistry::add<HyperboloidShellGeometry>();
    ItemRegistry::add<MultiGaussianExpansionGeometry>();
    ItemRegistry::add<GenGeometry>();
    ItemRegistry::add<UniformBoxGeometry>();
    ItemRegistry::add<ReadFitsGeometry>();
    ItemRegistry::add<ImportedGeometry>();
    ItemRegistry::add<ParticleGeometry>();
    ItemRegistry::add<MeshGeometry>();
    ItemRegistry::add<AdaptiveMeshGeometry>();
    ItemRegistry::add<VoronoiMeshGeometry>();

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

    // spatial grids
    ItemRegistry::add<SpatialGrid>();
    ItemRegistry::add<SphereSpatialGrid>();
    ItemRegistry::add<Sphere1DSpatialGrid>();
    ItemRegistry::add<Sphere2DSpatialGrid>();
    ItemRegistry::add<CylinderSpatialGrid>();
    ItemRegistry::add<Cylinder2DSpatialGrid>();
    ItemRegistry::add<BoxSpatialGrid>();
    ItemRegistry::add<CartesianSpatialGrid>();

    // one-dimensional meshes for spatial grids
    ItemRegistry::add<Mesh>();
    ItemRegistry::add<MoveableMesh>();
    ItemRegistry::add<AnchoredMesh>();
    ItemRegistry::add<LinMesh>();
    ItemRegistry::add<PowMesh>();
    ItemRegistry::add<SymPowMesh>();
    ItemRegistry::add<LogMesh>();
    ItemRegistry::add<FileMesh>();

    // medium system and media
    ItemRegistry::add<MediumSystem>();
    ItemRegistry::add<Medium>();
    ItemRegistry::add<GeometricMedium>();

    // material normalizations
    ItemRegistry::add<MaterialNormalization>();
    ItemRegistry::add<MassMaterialNormalization>();
    ItemRegistry::add<AxisMaterialNormalization>();
    ItemRegistry::add<OpticalDepthMaterialNormalization>();
    ItemRegistry::add<MassColumnMaterialNormalization>();
    ItemRegistry::add<NumberColumnMaterialNormalization>();

    // material mixes
    ItemRegistry::add<MaterialMix>();
    ItemRegistry::add<SingleGrainDustMix>();
    ItemRegistry::add<MeanTrustBenchmarkDustMix>();
    ItemRegistry::add<ElectronMix>();

    // wavelength grids
    ItemRegistry::add<WavelengthGrid>();
    ItemRegistry::add<DisjointWavelengthGrid>();
    ItemRegistry::add<LogWavelengthGrid>();
    ItemRegistry::add<NestedLogWavelengthGrid>();
    ItemRegistry::add<BandWavelengthGrid>();
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

