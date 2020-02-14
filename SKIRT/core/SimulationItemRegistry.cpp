/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SimulationItemRegistry.hpp"
#include "ItemRegistry.hpp"
#include "SkirtUnitDef.hpp"

// ---> add new items below in alphabetical order

#include "AdaptiveMeshGeometry.hpp"
#include "AdaptiveMeshMedium.hpp"
#include "AdaptiveMeshSource.hpp"
#include "AdaptiveMeshSpatialGrid.hpp"
#include "AllCellsLibrary.hpp"
#include "AllSkyInstrument.hpp"
#include "BandLuminosityNormalization.hpp"
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
#include "ConfigurableBandWavelengthGrid.hpp"
#include "ConfigurableDustMix.hpp"
#include "ConicalAngularDistribution.hpp"
#include "ConicalShellGeometry.hpp"
#include "CrystalEnstatiteGrainComposition.hpp"
#include "CrystalForsteriteGrainComposition.hpp"
#include "CubicSplineSmoothingKernel.hpp"
#include "CubicalBackgroundSource.hpp"
#include "Cylinder2DSpatialGrid.hpp"
#include "CylindricalClipGeometryDecorator.hpp"
#include "CylindricalVectorField.hpp"
#include "DefaultDustTemperatureCutsProbe.hpp"
#include "DefaultMagneticFieldCutsProbe.hpp"
#include "DefaultMediaDensityCutsProbe.hpp"
#include "DefaultMediumVelocityCutsProbe.hpp"
#include "DefaultRadiationFieldCutsProbe.hpp"
#include "DensityTreePolicy.hpp"
#include "DiscreteWavelengthDistribution.hpp"
#include "DraineGraphiteGrainComposition.hpp"
#include "DraineIonizedPAHGrainComposition.hpp"
#include "DraineLiDustMix.hpp"
#include "DraineNeutralPAHGrainComposition.hpp"
#include "DraineSilicateGrainComposition.hpp"
#include "DustAbsorptionPerCellProbe.hpp"
#include "DustEmGrainComposition.hpp"
#include "DustEmissionOptions.hpp"
#include "DustEmissionWavelengthGridProbe.hpp"
#include "DustEmissivityProbe.hpp"
#include "DustGrainPopulationsProbe.hpp"
#include "DustGrainSizeDistributionProbe.hpp"
#include "DustSelfAbsorptionOptions.hpp"
#include "DustTemperaturePerCellProbe.hpp"
#include "EinastoGeometry.hpp"
#include "ElectronMix.hpp"
#include "ExpDiskGeometry.hpp"
#include "ExtinctionOnlyOptions.hpp"
#include "ExtragalacticUnits.hpp"
#include "FSPSSED.hpp"
#include "FSPSSEDFamily.hpp"
#include "FieldStrengthCellLibrary.hpp"
#include "FileBand.hpp"
#include "FileMesh.hpp"
#include "FileSED.hpp"
#include "FileSSPSEDFamily.hpp"
#include "FileTreeSpatialGrid.hpp"
#include "FileWavelengthDistribution.hpp"
#include "FileWavelengthGrid.hpp"
#include "FlatUniverseCosmology.hpp"
#include "FrameInstrument.hpp"
#include "FullInstrument.hpp"
#include "GammaGeometry.hpp"
#include "GasMix.hpp"
#include "GaussianGeometry.hpp"
#include "GeometricMedium.hpp"
#include "GeometricSource.hpp"
#include "GrainPopulation.hpp"
#include "HammerAitoffProjection.hpp"
#include "HyperboloidGeometry.hpp"
#include "HyperboloidShellGeometry.hpp"
#include "InstrumentSystem.hpp"
#include "InstrumentWavelengthGridProbe.hpp"
#include "IntegratedLuminosityNormalization.hpp"
#include "IsotropicAngularDistribution.hpp"
#include "LaserAngularDistribution.hpp"
#include "LaunchedPacketsProbe.hpp"
#include "LinMesh.hpp"
#include "LinWavelengthDistribution.hpp"
#include "LinWavelengthGrid.hpp"
#include "LinearDustTemperatureCutProbe.hpp"
#include "ListBand.hpp"
#include "ListMesh.hpp"
#include "ListSED.hpp"
#include "ListWavelengthDistribution.hpp"
#include "ListWavelengthGrid.hpp"
#include "LocalUniverseCosmology.hpp"
#include "LogMesh.hpp"
#include "LogWavelengthDistribution.hpp"
#include "LogWavelengthGrid.hpp"
#include "LuminosityProbe.hpp"
#include "LyaDoublePeakedSED.hpp"
#include "LyaDoublePeakedSEDFamily.hpp"
#include "LyaGaussianSED.hpp"
#include "LyaGaussianSEDFamily.hpp"
#include "LyaSEDDecorator.hpp"
#include "LyaSEDFamilyDecorator.hpp"
#include "MRNDustMix.hpp"
#include "MagneticFieldPerCellProbe.hpp"
#include "MappingsSED.hpp"
#include "MappingsSEDFamily.hpp"
#include "MarastonSED.hpp"
#include "MarastonSEDFamily.hpp"
#include "MassColumnMaterialNormalization.hpp"
#include "MassMaterialNormalization.hpp"
#include "MeanFileDustMix.hpp"
#include "MeanInterstellarDustMix.hpp"
#include "MeanIvezicBenchmarkDustMix.hpp"
#include "MeanListDustMix.hpp"
#include "MeanPascucciBenchmarkDustMix.hpp"
#include "MeanPinteBenchmarkDustMix.hpp"
#include "MeanTrustBenchmarkDustMix.hpp"
#include "MediumSystem.hpp"
#include "MediumVelocityPerCellProbe.hpp"
#include "MeridionalDustTemperatureCutProbe.hpp"
#include "MieSilicateGrainComposition.hpp"
#include "MinSilicateGrainComposition.hpp"
#include "ModifiedLogNormalGrainSizeDistribution.hpp"
#include "ModifiedPowerLawGrainSizeDistribution.hpp"
#include "MollweideProjection.hpp"
#include "MonteCarloSimulation.hpp"
#include "MultiGaussianExpansionGeometry.hpp"
#include "NestedLogWavelengthGrid.hpp"
#include "NetzerAngularDistribution.hpp"
#include "NoPolarizationProfile.hpp"
#include "NumberColumnMaterialNormalization.hpp"
#include "OffsetGeometryDecorator.hpp"
#include "OffsetVectorFieldDecorator.hpp"
#include "OpticalDepthMapProbe.hpp"
#include "OpticalDepthMaterialNormalization.hpp"
#include "OpticalMaterialPropertiesProbe.hpp"
#include "ParaboloidGeometry.hpp"
#include "ParaboloidShellGeometry.hpp"
#include "ParticleGeometry.hpp"
#include "ParticleMedium.hpp"
#include "ParticleSource.hpp"
#include "PerspectiveInstrument.hpp"
#include "PhotonPacketOptions.hpp"
#include "PlanarDustTemperatureCutsProbe.hpp"
#include "PlanarMagneticFieldCutsProbe.hpp"
#include "PlanarMediaDensityCutsProbe.hpp"
#include "PlanarMediumVelocityCutsProbe.hpp"
#include "PlanarRadiationFieldCutsProbe.hpp"
#include "PlummerGeometry.hpp"
#include "PointSource.hpp"
#include "PolarizedGraphiteGrainComposition.hpp"
#include "PolarizedSilicateGrainComposition.hpp"
#include "PolicyTreeSpatialGrid.hpp"
#include "PowMesh.hpp"
#include "PowerLawGrainSizeDistribution.hpp"
#include "PredefinedBandWavelengthGrid.hpp"
#include "ProbeSystem.hpp"
#include "PseudoSersicGeometry.hpp"
#include "QuasarSED.hpp"
#include "RadialVectorField.hpp"
#include "RadiationFieldPerCellProbe.hpp"
#include "RadiationFieldWavelengthGridProbe.hpp"
#include "Random.hpp"
#include "ReadFitsGeometry.hpp"
#include "RingGeometry.hpp"
#include "RotateGeometryDecorator.hpp"
#include "RotateVectorFieldDecorator.hpp"
#include "SEDInstrument.hpp"
#include "SIUnits.hpp"
#include "ScaledGaussianSmoothingKernel.hpp"
#include "SelectDustMixFamily.hpp"
#include "SersicGeometry.hpp"
#include "ShellGeometry.hpp"
#include "SineSquarePolarizationProfile.hpp"
#include "SingleGrainSizeDistribution.hpp"
#include "SiteListTreePolicy.hpp"
#include "SourceSystem.hpp"
#include "SpatialCellPropertiesProbe.hpp"
#include "SpatialGrid.hpp"
#include "SpatialGridConvergenceProbe.hpp"
#include "SpatialGridPlotProbe.hpp"
#include "SpatialGridSourceDensityProbe.hpp"
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
#include "TemperatureWavelengthCellLibrary.hpp"
#include "ThemisDustMix.hpp"
#include "TorusGeometry.hpp"
#include "TreePolicy.hpp"
#include "TreeSpatialGrid.hpp"
#include "TreeSpatialGridTopologyProbe.hpp"
#include "TriaxialGeometryDecorator.hpp"
#include "TrustBenchmarkDustMix.hpp"
#include "TrustGraphiteGrainComposition.hpp"
#include "TrustNeutralPAHGrainComposition.hpp"
#include "TrustSilicateGrainComposition.hpp"
#include "UnidirectionalVectorField.hpp"
#include "UniformBoxGeometry.hpp"
#include "UniformSmoothingKernel.hpp"
#include "VoronoiMeshGeometry.hpp"
#include "VoronoiMeshMedium.hpp"
#include "VoronoiMeshSource.hpp"
#include "VoronoiMeshSpatialGrid.hpp"
#include "WeingartnerDraineDustMix.hpp"
#include "ZubkoDustMix.hpp"
#include "ZubkoGraphiteGrainSizeDistribution.hpp"
#include "ZubkoPAHGrainSizeDistribution.hpp"
#include "ZubkoSilicateGrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

SimulationItemRegistry::SimulationItemRegistry(string version, string format)
{
    // start a new schema
    ItemRegistry::beginSchema("SKIRT", "a SKIRT parameter file", version, "ski", "skirt-simulation-hierarchy",
                              "MonteCarloSimulation", format, "http://www.skirt.ugent.be/skirt9");

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
    ItemRegistry::add<Cosmology>();
    ItemRegistry::add<LocalUniverseCosmology>();
    ItemRegistry::add<FlatUniverseCosmology>();

    // source system and sources
    ItemRegistry::add<SourceSystem>();
    ItemRegistry::add<Source>();
    ItemRegistry::add<NormalizedSource>();
    ItemRegistry::add<SpecialtySource>();
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
    ItemRegistry::add<FSPSSED>();
    ItemRegistry::add<MappingsSED>();
    ItemRegistry::add<TabulatedSED>();
    ItemRegistry::add<FileSED>();
    ItemRegistry::add<ListSED>();
    ItemRegistry::add<LyaGaussianSED>();
    ItemRegistry::add<LyaDoublePeakedSED>();
    ItemRegistry::add<LyaSEDDecorator>();

    // SED families
    ItemRegistry::add<SEDFamily>();
    ItemRegistry::add<BlackBodySEDFamily>();
    ItemRegistry::add<CastelliKuruczSEDFamily>();
    ItemRegistry::add<BruzualCharlotSEDFamily>();
    ItemRegistry::add<MarastonSEDFamily>();
    ItemRegistry::add<Starburst99SEDFamily>();
    ItemRegistry::add<FSPSSEDFamily>();
    ItemRegistry::add<FileSSPSEDFamily>();
    ItemRegistry::add<MappingsSEDFamily>();
    ItemRegistry::add<LyaGaussianSEDFamily>();
    ItemRegistry::add<LyaDoublePeakedSEDFamily>();
    ItemRegistry::add<LyaSEDFamilyDecorator>();

    // wavelength distributions
    ItemRegistry::add<WavelengthDistribution>();
    ItemRegistry::add<RangeWavelengthDistribution>();
    ItemRegistry::add<LinWavelengthDistribution>();
    ItemRegistry::add<LogWavelengthDistribution>();
    ItemRegistry::add<TabulatedWavelengthDistribution>();
    ItemRegistry::add<FileWavelengthDistribution>();
    ItemRegistry::add<ListWavelengthDistribution>();
    ItemRegistry::add<DiscreteWavelengthDistribution>();

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

    // vector fields
    ItemRegistry::add<VectorField>();
    ItemRegistry::add<RadialVectorField>();
    ItemRegistry::add<CylindricalVectorField>();
    ItemRegistry::add<UnidirectionalVectorField>();
    ItemRegistry::add<OffsetVectorFieldDecorator>();
    ItemRegistry::add<RotateVectorFieldDecorator>();

    // spatial grids
    ItemRegistry::add<SpatialGrid>();
    ItemRegistry::add<SphereSpatialGrid>();
    ItemRegistry::add<Sphere1DSpatialGrid>();
    ItemRegistry::add<Sphere2DSpatialGrid>();
    ItemRegistry::add<CylinderSpatialGrid>();
    ItemRegistry::add<Cylinder2DSpatialGrid>();
    ItemRegistry::add<BoxSpatialGrid>();
    ItemRegistry::add<CartesianSpatialGrid>();
    ItemRegistry::add<TreeSpatialGrid>();
    ItemRegistry::add<PolicyTreeSpatialGrid>();
    ItemRegistry::add<FileTreeSpatialGrid>();
    ItemRegistry::add<AdaptiveMeshSpatialGrid>();
    ItemRegistry::add<VoronoiMeshSpatialGrid>();

    // spatial grid policies
    ItemRegistry::add<TreePolicy>();
    ItemRegistry::add<DensityTreePolicy>();
    ItemRegistry::add<SiteListTreePolicy>();

    // one-dimensional meshes for spatial grids
    ItemRegistry::add<Mesh>();
    ItemRegistry::add<MoveableMesh>();
    ItemRegistry::add<AnchoredMesh>();
    ItemRegistry::add<LinMesh>();
    ItemRegistry::add<PowMesh>();
    ItemRegistry::add<SymPowMesh>();
    ItemRegistry::add<LogMesh>();
    ItemRegistry::add<TabulatedMesh>();
    ItemRegistry::add<FileMesh>();
    ItemRegistry::add<ListMesh>();

    // medium system and media
    ItemRegistry::add<MediumSystem>();
    ItemRegistry::add<Medium>();
    ItemRegistry::add<GeometricMedium>();
    ItemRegistry::add<ImportedMedium>();
    ItemRegistry::add<ParticleMedium>();
    ItemRegistry::add<MeshMedium>();
    ItemRegistry::add<AdaptiveMeshMedium>();
    ItemRegistry::add<VoronoiMeshMedium>();

    // medium system options
    ItemRegistry::add<PhotonPacketOptions>();
    ItemRegistry::add<ExtinctionOnlyOptions>();
    ItemRegistry::add<DustEmissionOptions>();
    ItemRegistry::add<DustSelfAbsorptionOptions>();

    // material normalizations
    ItemRegistry::add<MaterialNormalization>();
    ItemRegistry::add<MassMaterialNormalization>();
    ItemRegistry::add<AxisMaterialNormalization>();
    ItemRegistry::add<OpticalDepthMaterialNormalization>();
    ItemRegistry::add<MassColumnMaterialNormalization>();
    ItemRegistry::add<NumberColumnMaterialNormalization>();

    // material mixes
    ItemRegistry::add<MaterialMix>();
    ItemRegistry::add<DustMix>();

    ItemRegistry::add<SingleGrainDustMix>();
    ItemRegistry::add<MeanInterstellarDustMix>();
    ItemRegistry::add<MeanTrustBenchmarkDustMix>();
    ItemRegistry::add<MeanPinteBenchmarkDustMix>();
    ItemRegistry::add<MeanPascucciBenchmarkDustMix>();
    ItemRegistry::add<MeanIvezicBenchmarkDustMix>();
    ItemRegistry::add<TabulatedDustMix>();
    ItemRegistry::add<MeanFileDustMix>();
    ItemRegistry::add<MeanListDustMix>();

    ItemRegistry::add<MultiGrainDustMix>();
    ItemRegistry::add<ThemisDustMix>();
    ItemRegistry::add<DraineLiDustMix>();
    ItemRegistry::add<ZubkoDustMix>();
    ItemRegistry::add<WeingartnerDraineDustMix>();
    ItemRegistry::add<MRNDustMix>();
    ItemRegistry::add<TrustBenchmarkDustMix>();
    ItemRegistry::add<ConfigurableDustMix>();

    ItemRegistry::add<ElectronMix>();

    ItemRegistry::add<GasMix>();

    // material mix families
    ItemRegistry::add<MaterialMixFamily>();
    ItemRegistry::add<SelectDustMixFamily>();

    // grain population
    ItemRegistry::add<GrainPopulation>();

    // grain size distributions
    ItemRegistry::add<GrainSizeDistribution>();
    ItemRegistry::add<RangeGrainSizeDistribution>();
    ItemRegistry::add<PowerLawGrainSizeDistribution>();
    ItemRegistry::add<ModifiedPowerLawGrainSizeDistribution>();
    ItemRegistry::add<LogNormalGrainSizeDistribution>();
    ItemRegistry::add<ModifiedLogNormalGrainSizeDistribution>();
    ItemRegistry::add<SingleGrainSizeDistribution>();
    ItemRegistry::add<ZubkoSilicateGrainSizeDistribution>();
    ItemRegistry::add<ZubkoGraphiteGrainSizeDistribution>();
    ItemRegistry::add<ZubkoPAHGrainSizeDistribution>();

    // grain compositions
    ItemRegistry::add<GrainComposition>();
    ItemRegistry::add<DraineSilicateGrainComposition>();
    ItemRegistry::add<DraineGraphiteGrainComposition>();
    ItemRegistry::add<DraineNeutralPAHGrainComposition>();
    ItemRegistry::add<DraineIonizedPAHGrainComposition>();
    ItemRegistry::add<MieSilicateGrainComposition>();
    ItemRegistry::add<MinSilicateGrainComposition>();
    ItemRegistry::add<PolarizedSilicateGrainComposition>();
    ItemRegistry::add<PolarizedGraphiteGrainComposition>();
    ItemRegistry::add<CrystalEnstatiteGrainComposition>();
    ItemRegistry::add<CrystalForsteriteGrainComposition>();
    ItemRegistry::add<DustEmGrainComposition>();
    ItemRegistry::add<TrustSilicateGrainComposition>();
    ItemRegistry::add<TrustGraphiteGrainComposition>();
    ItemRegistry::add<TrustNeutralPAHGrainComposition>();

    // spatial cell libraries
    ItemRegistry::add<SpatialCellLibrary>();
    ItemRegistry::add<AllCellsLibrary>();
    ItemRegistry::add<FieldStrengthCellLibrary>();
    ItemRegistry::add<TemperatureWavelengthCellLibrary>();

    // wavelength grids
    ItemRegistry::add<WavelengthGrid>();
    ItemRegistry::add<DisjointWavelengthGrid>();
    ItemRegistry::add<LogWavelengthGrid>();
    ItemRegistry::add<NestedLogWavelengthGrid>();
    ItemRegistry::add<LinWavelengthGrid>();
    ItemRegistry::add<BandWavelengthGrid>();
    ItemRegistry::add<PredefinedBandWavelengthGrid>();
    ItemRegistry::add<ConfigurableBandWavelengthGrid>();
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
    ItemRegistry::add<AbstractWavelengthProbe>();
    ItemRegistry::add<AbstractWavelengthGridProbe>();
    ItemRegistry::add<AbstractPlanarCutsProbe>();
    ItemRegistry::add<InstrumentWavelengthGridProbe>();
    ItemRegistry::add<LuminosityProbe>();
    ItemRegistry::add<LaunchedPacketsProbe>();
    ItemRegistry::add<SpatialGridPlotProbe>();
    ItemRegistry::add<SpatialGridConvergenceProbe>();
    ItemRegistry::add<TreeSpatialGridTopologyProbe>();
    ItemRegistry::add<DefaultMediaDensityCutsProbe>();
    ItemRegistry::add<PlanarMediaDensityCutsProbe>();
    ItemRegistry::add<OpticalDepthMapProbe>();
    ItemRegistry::add<SpatialCellPropertiesProbe>();
    ItemRegistry::add<SpatialGridSourceDensityProbe>();
    ItemRegistry::add<OpticalMaterialPropertiesProbe>();
    ItemRegistry::add<DustGrainPopulationsProbe>();
    ItemRegistry::add<DustGrainSizeDistributionProbe>();

    ItemRegistry::add<DefaultMediumVelocityCutsProbe>();
    ItemRegistry::add<PlanarMediumVelocityCutsProbe>();
    ItemRegistry::add<MediumVelocityPerCellProbe>();
    ItemRegistry::add<DefaultMagneticFieldCutsProbe>();
    ItemRegistry::add<PlanarMagneticFieldCutsProbe>();
    ItemRegistry::add<MagneticFieldPerCellProbe>();

    ItemRegistry::add<DefaultRadiationFieldCutsProbe>();
    ItemRegistry::add<PlanarRadiationFieldCutsProbe>();
    ItemRegistry::add<RadiationFieldPerCellProbe>();
    ItemRegistry::add<RadiationFieldWavelengthGridProbe>();
    ItemRegistry::add<DefaultDustTemperatureCutsProbe>();
    ItemRegistry::add<PlanarDustTemperatureCutsProbe>();
    ItemRegistry::add<DustTemperaturePerCellProbe>();
    ItemRegistry::add<LinearDustTemperatureCutProbe>();
    ItemRegistry::add<MeridionalDustTemperatureCutProbe>();
    ItemRegistry::add<DustAbsorptionPerCellProbe>();
    ItemRegistry::add<DustEmissivityProbe>();
    ItemRegistry::add<DustEmissionWavelengthGridProbe>();

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
