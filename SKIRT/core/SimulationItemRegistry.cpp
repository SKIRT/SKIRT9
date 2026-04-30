/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SimulationItemRegistry.hpp"
#include "ItemRegistry.hpp"
#include "SkirtUnitDef.hpp"

// ---> add new items below in alphabetical order

#include "AbsorptionOnlyMaterialMixDecorator.hpp"
#include "AdaptiveMeshGeometry.hpp"
#include "AdaptiveMeshMedium.hpp"
#include "AdaptiveMeshSource.hpp"
#include "AdaptiveMeshSpatialGrid.hpp"
#include "AllCellsLibrary.hpp"
#include "AllSkyInstrument.hpp"
#include "AllSkyProjectionForm.hpp"
#include "AnnulusGeometry.hpp"
#include "AtPositionsForm.hpp"
#include "AxPowerLawRedistributeGeometryDecorator.hpp"
#include "BandLuminosityNormalization.hpp"
#include "BegemannPorousAluminaGrainComposition.hpp"
#include "BlackBodySED.hpp"
#include "BlackBodySEDFamily.hpp"
#include "BoxClipGeometryDecorator.hpp"
#include "BpassSED.hpp"
#include "BpassSEDFamily.hpp"
#include "BroadBand.hpp"
#include "BrokenExpDiskGeometry.hpp"
#include "BruzualCharlotSED.hpp"
#include "BruzualCharlotSEDFamily.hpp"
#include "CartesianSpatialGrid.hpp"
#include "CastelliKuruczSED.hpp"
#include "CastelliKuruczSEDFamily.hpp"
#include "CellGeometry.hpp"
#include "CellMedium.hpp"
#include "CellSource.hpp"
#include "ClearDensityRecipe.hpp"
#include "ClumpyGeometryDecorator.hpp"
#include "CombineGeometryDecorator.hpp"
#include "CompositeWavelengthGrid.hpp"
#include "ConfigurableBandWavelengthGrid.hpp"
#include "ConfigurableDustMix.hpp"
#include "ConicalAngularDistribution.hpp"
#include "ConicalShellGeometry.hpp"
#include "ConvergenceCutsProbe.hpp"
#include "ConvergenceInfoProbe.hpp"
#include "CrystalEnstatiteGrainComposition.hpp"
#include "CrystalForsteriteGrainComposition.hpp"
#include "CubicSplineSmoothingKernel.hpp"
#include "CubicalBackgroundSource.hpp"
#include "CustomStateProbe.hpp"
#include "Cylinder2DSpatialGrid.hpp"
#include "Cylinder3DSpatialGrid.hpp"
#include "CylindricalCellGeometry.hpp"
#include "CylindricalCellMedium.hpp"
#include "CylindricalCellSource.hpp"
#include "CylindricalClipGeometryDecorator.hpp"
#include "CylindricalVectorField.hpp"
#include "DefaultCutsForm.hpp"
#include "DefaultWavelengthDistribution.hpp"
#include "DensityProbe.hpp"
#include "DensityTreePolicy.hpp"
#include "DiscreteWavelengthDistribution.hpp"
#include "DonutGeometry.hpp"
#include "DorschnerOlivineGrainComposition.hpp"
#include "DraineGraphiteGrainComposition.hpp"
#include "DraineIonizedPAHGrainComposition.hpp"
#include "DraineLiDustMix.hpp"
#include "DraineNeutralPAHGrainComposition.hpp"
#include "DraineSilicateGrainComposition.hpp"
#include "DustAbsorptionPerCellProbe.hpp"
#include "DustEmGrainComposition.hpp"
#include "DustEmissionWavelengthGridProbe.hpp"
#include "DustEmissivityProbe.hpp"
#include "DustGrainPopulationsProbe.hpp"
#include "DustGrainSizeDistributionProbe.hpp"
#include "EinastoGeometry.hpp"
#include "ElectronMix.hpp"
#include "ExpDiskGeometry.hpp"
#include "ExtragalacticUnits.hpp"
#include "FSPSSED.hpp"
#include "FSPSSEDFamily.hpp"
#include "FieldStrengthCellLibrary.hpp"
#include "FileBand.hpp"
#include "FileBorderWavelengthGrid.hpp"
#include "FileGrainSizeDistribution.hpp"
#include "FileIndexedSEDFamily.hpp"
#include "FileLineSED.hpp"
#include "FileMesh.hpp"
#include "FilePolarizedPointSource.hpp"
#include "FileSED.hpp"
#include "FileSSPSEDFamily.hpp"
#include "FileTreeSpatialGrid.hpp"
#include "FileWavelengthDistribution.hpp"
#include "FileWavelengthGrid.hpp"
#include "FlatUniverseCosmology.hpp"
#include "FragmentDustMixDecorator.hpp"
#include "FrameInstrument.hpp"
#include "FullInstrument.hpp"
#include "GammaGeometry.hpp"
#include "GaussianGeometry.hpp"
#include "GeometricMedium.hpp"
#include "GeometricSource.hpp"
#include "GrainPopulation.hpp"
#include "HEALPixSkyInstrument.hpp"
#include "HammerAitoffProjection.hpp"
#include "HirashitaLogNormalGrainSizeDistribution.hpp"
#include "HofmeisterPericlaseGrainComposition.hpp"
#include "HollowRadialVectorField.hpp"
#include "HubbleRadialVectorField.hpp"
#include "HyperboloidGeometry.hpp"
#include "HyperboloidShellGeometry.hpp"
#include "ImportedMediumDensityProbe.hpp"
#include "ImportedMediumMetallicityProbe.hpp"
#include "ImportedMediumTemperatureProbe.hpp"
#include "ImportedMediumVelocityProbe.hpp"
#include "ImportedSourceAgeProbe.hpp"
#include "ImportedSourceDensityProbe.hpp"
#include "ImportedSourceLuminosityProbe.hpp"
#include "ImportedSourceMetallicityProbe.hpp"
#include "ImportedSourceVelocityProbe.hpp"
#include "InstrumentSystem.hpp"
#include "InstrumentWavelengthGridProbe.hpp"
#include "IntegratedLuminosityNormalization.hpp"
#include "IsotropicAngularDistribution.hpp"
#include "LaserAngularDistribution.hpp"
#include "LaunchedPacketsProbe.hpp"
#include "LinBorderWavelengthGrid.hpp"
#include "LinMesh.hpp"
#include "LinWavelengthDistribution.hpp"
#include "LinWavelengthGrid.hpp"
#include "LineLuminosityNormalization.hpp"
#include "LinearCutForm.hpp"
#include "LinearDustDestructionRecipe.hpp"
#include "ListBand.hpp"
#include "ListBorderWavelengthGrid.hpp"
#include "ListGrainSizeDistribution.hpp"
#include "ListLineSED.hpp"
#include "ListMesh.hpp"
#include "ListSED.hpp"
#include "ListWavelengthDistribution.hpp"
#include "ListWavelengthGrid.hpp"
#include "LocalUniverseCosmology.hpp"
#include "LogBorderWavelengthGrid.hpp"
#include "LogMesh.hpp"
#include "LogWavelengthDistribution.hpp"
#include "LogWavelengthGrid.hpp"
#include "LuminosityProbe.hpp"
#include "LyaDoublePeakedSED.hpp"
#include "LyaDoublePeakedSEDFamily.hpp"
#include "LyaGaussianSED.hpp"
#include "LyaGaussianSEDFamily.hpp"
#include "LyaNeutralHydrogenGasMix.hpp"
#include "LyaSEDDecorator.hpp"
#include "LyaSEDFamilyDecorator.hpp"
#include "MRNDustMix.hpp"
#include "MagneticFieldProbe.hpp"
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
#include "MeridionalCutForm.hpp"
#include "MetallicityProbe.hpp"
#include "MieSilicateGrainComposition.hpp"
#include "MinSilicateGrainComposition.hpp"
#include "ModifiedLogNormalGrainSizeDistribution.hpp"
#include "ModifiedPowerLawGrainSizeDistribution.hpp"
#include "MollweideProjection.hpp"
#include "MonteCarloSimulation.hpp"
#include "MultiGaussianExpansionGeometry.hpp"
#include "NestedDensityTreePolicy.hpp"
#include "NestedLogWavelengthGrid.hpp"
#include "NetzerAngularDistribution.hpp"
#include "NoPolarizationProfile.hpp"
#include "NonLTELineGasMix.hpp"
#include "NumberColumnMaterialNormalization.hpp"
#include "NumberMaterialNormalization.hpp"
#include "OffsetGeometryDecorator.hpp"
#include "OffsetVectorFieldDecorator.hpp"
#include "OpacityProbe.hpp"
#include "OpticalDepthMaterialNormalization.hpp"
#include "OpticalMaterialPropertiesProbe.hpp"
#include "ParaboloidGeometry.hpp"
#include "ParaboloidShellGeometry.hpp"
#include "ParallelProjectionForm.hpp"
#include "ParticleGeometry.hpp"
#include "ParticleMedium.hpp"
#include "ParticleSource.hpp"
#include "PerCellForm.hpp"
#include "PerspectiveInstrument.hpp"
#include "PlanarCutsForm.hpp"
#include "PlummerGeometry.hpp"
#include "PointSource.hpp"
#include "PolicyTreeSpatialGrid.hpp"
#include "PowMesh.hpp"
#include "PowerLawGrainSizeDistribution.hpp"
#include "PredefinedBandWavelengthGrid.hpp"
#include "ProbeSystem.hpp"
#include "PseudoSersicGeometry.hpp"
#include "QuarticSplineSmoothingKernel.hpp"
#include "QuasarSED.hpp"
#include "RadialVectorField.hpp"
#include "RadiationFieldProbe.hpp"
#include "RadiationFieldWavelengthGridProbe.hpp"
#include "Random.hpp"
#include "ReadFits3DGeometry.hpp"
#include "ReadFitsGeometry.hpp"
#include "RedistributeGeometryDecorator.hpp"
#include "ResolutionBorderWavelengthGrid.hpp"
#include "ResolutionWavelengthGrid.hpp"
#include "RingGeometry.hpp"
#include "RotateGeometryDecorator.hpp"
#include "RotateVectorFieldDecorator.hpp"
#include "SEDInstrument.hpp"
#include "SIUnits.hpp"
#include "ScaledGaussianSmoothingKernel.hpp"
#include "SecondaryDustLuminosityProbe.hpp"
#include "SecondaryLineLuminosityProbe.hpp"
#include "SelectDustMixFamily.hpp"
#include "SersicGeometry.hpp"
#include "ShellGeometry.hpp"
#include "SineSquarePolarizationProfile.hpp"
#include "SingleGrainSizeDistribution.hpp"
#include "SingleWavelengthSED.hpp"
#include "SiteListTreePolicy.hpp"
#include "SourceSystem.hpp"
#include "SpatialCellPropertiesProbe.hpp"
#include "SpatialGrid.hpp"
#include "SpatialGridPlotProbe.hpp"
#include "SpatialGridSourceDensityProbe.hpp"
#include "SpecificLuminosityNormalization.hpp"
#include "SphePowerLawRedistributeGeometryDecorator.hpp"
#include "Sphere1DSpatialGrid.hpp"
#include "Sphere2DSpatialGrid.hpp"
#include "Sphere3DSpatialGrid.hpp"
#include "SphericalBackgroundSource.hpp"
#include "SphericalCellGeometry.hpp"
#include "SphericalCellMedium.hpp"
#include "SphericalCellSource.hpp"
#include "SphericalClipGeometryDecorator.hpp"
#include "SpheroidalGeometryDecorator.hpp"
#include "SpheroidalGraphiteGrainComposition.hpp"
#include "SpheroidalSilicateGrainComposition.hpp"
#include "SpinFlipAbsorptionMix.hpp"
#include "SpinFlipHydrogenGasMix.hpp"
#include "SpinFlipSEDFamily.hpp"
#include "SpiralStructureGeometryDecorator.hpp"
#include "Starburst99SED.hpp"
#include "Starburst99SEDFamily.hpp"
#include "StellarSurfaceSource.hpp"
#include "StellarUnits.hpp"
#include "SunSED.hpp"
#include "SymCosMesh.hpp"
#include "SymLogMesh.hpp"
#include "SymPowMesh.hpp"
#include "TTauriDiskGeometry.hpp"
#include "TemperatureProbe.hpp"
#include "TemperatureWavelengthCellLibrary.hpp"
#include "TetraMeshSpatialGrid.hpp"
#include "ThemisDustMix.hpp"
#include "ToddlersSED.hpp"
#include "ToddlersSEDFamily.hpp"
#include "TorusGeometry.hpp"
#include "TreePolicy.hpp"
#include "TreeSpatialGrid.hpp"
#include "TreeSpatialGridTopologyProbe.hpp"
#include "TriaxialGeometryDecorator.hpp"
#include "TrivialGasMix.hpp"
#include "TrustBenchmarkDustMix.hpp"
#include "TrustGraphiteGrainComposition.hpp"
#include "TrustNeutralPAHGrainComposition.hpp"
#include "TrustSilicateGrainComposition.hpp"
#include "UnidirectionalVectorField.hpp"
#include "UniformBoxGeometry.hpp"
#include "UniformSmoothingKernel.hpp"
#include "VelocityProbe.hpp"
#include "VoronoiMeshGeometry.hpp"
#include "VoronoiMeshMedium.hpp"
#include "VoronoiMeshSource.hpp"
#include "VoronoiMeshSpatialGrid.hpp"
#include "WeingartnerDraineDustMix.hpp"
#include "XRayAtomicGasMix.hpp"
#include "XRayIonicGasMix.hpp"
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
    ItemRegistry::add<CellSource>();
    ItemRegistry::add<CylindricalCellSource>();
    ItemRegistry::add<SphericalCellSource>();
    ItemRegistry::add<MeshSource>();
    ItemRegistry::add<AdaptiveMeshSource>();
    ItemRegistry::add<VoronoiMeshSource>();
    ItemRegistry::add<CenteredSource>();
    ItemRegistry::add<StellarSurfaceSource>();
    ItemRegistry::add<CubicalBackgroundSource>();
    ItemRegistry::add<SphericalBackgroundSource>();
    ItemRegistry::add<FilePolarizedPointSource>();

    // luminosity normalizations
    ItemRegistry::add<LuminosityNormalization>();
    ItemRegistry::add<IntegratedLuminosityNormalization>();
    ItemRegistry::add<SpecificLuminosityNormalization>();
    ItemRegistry::add<BandLuminosityNormalization>();
    ItemRegistry::add<LineLuminosityNormalization>();

    // SEDs
    ItemRegistry::add<SED>();
    ItemRegistry::add<ContSED>();
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
    ItemRegistry::add<BpassSED>();
    ItemRegistry::add<MappingsSED>();
    ItemRegistry::add<ToddlersSED>();
    ItemRegistry::add<TabulatedSED>();
    ItemRegistry::add<FileSED>();
    ItemRegistry::add<ListSED>();
    ItemRegistry::add<LyaGaussianSED>();
    ItemRegistry::add<LyaDoublePeakedSED>();
    ItemRegistry::add<LyaSEDDecorator>();
    ItemRegistry::add<LineSED>();
    ItemRegistry::add<FileLineSED>();
    ItemRegistry::add<ListLineSED>();
    ItemRegistry::add<SingleWavelengthSED>();

    // SED families
    ItemRegistry::add<SEDFamily>();
    ItemRegistry::add<BlackBodySEDFamily>();
    ItemRegistry::add<CastelliKuruczSEDFamily>();
    ItemRegistry::add<BruzualCharlotSEDFamily>();
    ItemRegistry::add<MarastonSEDFamily>();
    ItemRegistry::add<Starburst99SEDFamily>();
    ItemRegistry::add<FSPSSEDFamily>();
    ItemRegistry::add<BpassSEDFamily>();
    ItemRegistry::add<FileSSPSEDFamily>();
    ItemRegistry::add<FileIndexedSEDFamily>();
    ItemRegistry::add<MappingsSEDFamily>();
    ItemRegistry::add<ToddlersSEDFamily>();
    ItemRegistry::add<SpinFlipSEDFamily>();
    ItemRegistry::add<LyaGaussianSEDFamily>();
    ItemRegistry::add<LyaDoublePeakedSEDFamily>();
    ItemRegistry::add<LyaSEDFamilyDecorator>();

    // wavelength distributions
    ItemRegistry::add<WavelengthDistribution>();
    ItemRegistry::add<DefaultWavelengthDistribution>();
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
    ItemRegistry::add<TabulatedBand>();
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
    ItemRegistry::add<DonutGeometry>();
    ItemRegistry::add<AnnulusGeometry>();
    ItemRegistry::add<TTauriDiskGeometry>();
    ItemRegistry::add<ConicalShellGeometry>();
    ItemRegistry::add<ParaboloidGeometry>();
    ItemRegistry::add<ParaboloidShellGeometry>();
    ItemRegistry::add<HyperboloidGeometry>();
    ItemRegistry::add<HyperboloidShellGeometry>();
    ItemRegistry::add<MultiGaussianExpansionGeometry>();
    ItemRegistry::add<GenGeometry>();
    ItemRegistry::add<UniformBoxGeometry>();
    ItemRegistry::add<ReadFits3DGeometry>();
    ItemRegistry::add<ReadFitsGeometry>();
    ItemRegistry::add<ImportedGeometry>();
    ItemRegistry::add<ParticleGeometry>();
    ItemRegistry::add<CellGeometry>();
    ItemRegistry::add<CylindricalCellGeometry>();
    ItemRegistry::add<SphericalCellGeometry>();
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
    ItemRegistry::add<RedistributeGeometryDecorator>();
    ItemRegistry::add<AxPowerLawRedistributeGeometryDecorator>();
    ItemRegistry::add<SphePowerLawRedistributeGeometryDecorator>();

    // smoothing kernels
    ItemRegistry::add<SmoothingKernel>();
    ItemRegistry::add<CubicSplineSmoothingKernel>();
    ItemRegistry::add<QuarticSplineSmoothingKernel>();
    ItemRegistry::add<ScaledGaussianSmoothingKernel>();
    ItemRegistry::add<UniformSmoothingKernel>();

    // vector fields
    ItemRegistry::add<VectorField>();
    ItemRegistry::add<RadialVectorField>();
    ItemRegistry::add<HollowRadialVectorField>();
    ItemRegistry::add<HubbleRadialVectorField>();
    ItemRegistry::add<CylindricalVectorField>();
    ItemRegistry::add<UnidirectionalVectorField>();
    ItemRegistry::add<OffsetVectorFieldDecorator>();
    ItemRegistry::add<RotateVectorFieldDecorator>();

    // spatial grids
    ItemRegistry::add<SpatialGrid>();
    ItemRegistry::add<SphereSpatialGrid>();
    ItemRegistry::add<CylinderSpatialGrid>();
    ItemRegistry::add<Sphere1DSpatialGrid>();
    ItemRegistry::add<Sphere2DSpatialGrid>();
    ItemRegistry::add<Cylinder2DSpatialGrid>();
    ItemRegistry::add<Sphere3DSpatialGrid>();
    ItemRegistry::add<Cylinder3DSpatialGrid>();
    ItemRegistry::add<BoxSpatialGrid>();
    ItemRegistry::add<CartesianSpatialGrid>();
    ItemRegistry::add<TreeSpatialGrid>();
    ItemRegistry::add<PolicyTreeSpatialGrid>();
    ItemRegistry::add<FileTreeSpatialGrid>();
    ItemRegistry::add<AdaptiveMeshSpatialGrid>();
    ItemRegistry::add<VoronoiMeshSpatialGrid>();
    ItemRegistry::add<TetraMeshSpatialGrid>();

    // spatial grid policies
    ItemRegistry::add<TreePolicy>();
    ItemRegistry::add<DensityTreePolicy>();
    ItemRegistry::add<NestedDensityTreePolicy>();
    ItemRegistry::add<SiteListTreePolicy>();

    // one-dimensional meshes for spatial grids
    ItemRegistry::add<Mesh>();
    ItemRegistry::add<LinMesh>();
    ItemRegistry::add<PowMesh>();
    ItemRegistry::add<SymPowMesh>();
    ItemRegistry::add<LogMesh>();
    ItemRegistry::add<SymLogMesh>();
    ItemRegistry::add<SymCosMesh>();
    ItemRegistry::add<TabulatedMesh>();
    ItemRegistry::add<FileMesh>();
    ItemRegistry::add<ListMesh>();

    // medium system and media
    ItemRegistry::add<MediumSystem>();
    ItemRegistry::add<Medium>();
    ItemRegistry::add<GeometricMedium>();
    ItemRegistry::add<ImportedMedium>();
    ItemRegistry::add<ParticleMedium>();
    ItemRegistry::add<CellMedium>();
    ItemRegistry::add<CylindricalCellMedium>();
    ItemRegistry::add<SphericalCellMedium>();
    ItemRegistry::add<MeshMedium>();
    ItemRegistry::add<AdaptiveMeshMedium>();
    ItemRegistry::add<VoronoiMeshMedium>();

    // medium system options
    ItemRegistry::add<PhotonPacketOptions>();
    ItemRegistry::add<LyaOptions>();
    ItemRegistry::add<DynamicStateOptions>();
    ItemRegistry::add<RadiationFieldOptions>();
    ItemRegistry::add<SecondaryEmissionOptions>();
    ItemRegistry::add<IterationOptions>();
    ItemRegistry::add<DustEmissionOptions>();
    ItemRegistry::add<SamplingOptions>();

    // material normalizations
    ItemRegistry::add<MaterialNormalization>();
    ItemRegistry::add<MassMaterialNormalization>();
    ItemRegistry::add<NumberMaterialNormalization>();
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

    ItemRegistry::add<FragmentDustMixDecorator>();

    ItemRegistry::add<ElectronMix>();
    ItemRegistry::add<SpinFlipAbsorptionMix>();
    ItemRegistry::add<SpinFlipHydrogenGasMix>();
    ItemRegistry::add<XRayAtomicGasMix>();
    ItemRegistry::add<XRayIonicGasMix>();
    ItemRegistry::add<EmittingGasMix>();
    ItemRegistry::add<NonLTELineGasMix>();
    ItemRegistry::add<LyaNeutralHydrogenGasMix>();
    ItemRegistry::add<TrivialGasMix>();

    ItemRegistry::add<AbsorptionOnlyMaterialMixDecorator>();

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
    ItemRegistry::add<HirashitaLogNormalGrainSizeDistribution>();
    ItemRegistry::add<FileGrainSizeDistribution>();
    ItemRegistry::add<ListGrainSizeDistribution>();

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
    ItemRegistry::add<BegemannPorousAluminaGrainComposition>();
    ItemRegistry::add<HofmeisterPericlaseGrainComposition>();
    ItemRegistry::add<DorschnerOlivineGrainComposition>();
    ItemRegistry::add<TrustSilicateGrainComposition>();
    ItemRegistry::add<TrustGraphiteGrainComposition>();
    ItemRegistry::add<TrustNeutralPAHGrainComposition>();
    ItemRegistry::add<SpheroidalSilicateGrainComposition>();
    ItemRegistry::add<SpheroidalGraphiteGrainComposition>();

    // spatial cell libraries
    ItemRegistry::add<SpatialCellLibrary>();
    ItemRegistry::add<AllCellsLibrary>();
    ItemRegistry::add<FieldStrengthCellLibrary>();
    ItemRegistry::add<TemperatureWavelengthCellLibrary>();

    // dynamic medium state recipes
    ItemRegistry::add<DynamicStateRecipe>();
    ItemRegistry::add<ClearDensityRecipe>();
    ItemRegistry::add<DustDestructionRecipe>();
    ItemRegistry::add<LinearDustDestructionRecipe>();

    // wavelength grids
    ItemRegistry::add<WavelengthGrid>();
    ItemRegistry::add<DisjointWavelengthGrid>();
    ItemRegistry::add<LogWavelengthGrid>();
    ItemRegistry::add<NestedLogWavelengthGrid>();
    ItemRegistry::add<LinWavelengthGrid>();
    ItemRegistry::add<ResolutionWavelengthGrid>();
    ItemRegistry::add<FileWavelengthGrid>();
    ItemRegistry::add<ListWavelengthGrid>();
    ItemRegistry::add<LogBorderWavelengthGrid>();
    ItemRegistry::add<LinBorderWavelengthGrid>();
    ItemRegistry::add<ResolutionBorderWavelengthGrid>();
    ItemRegistry::add<FileBorderWavelengthGrid>();
    ItemRegistry::add<ListBorderWavelengthGrid>();
    ItemRegistry::add<CompositeWavelengthGrid>();
    ItemRegistry::add<BandWavelengthGrid>();
    ItemRegistry::add<PredefinedBandWavelengthGrid>();
    ItemRegistry::add<ConfigurableBandWavelengthGrid>();

    // instrument system and instruments
    ItemRegistry::add<InstrumentSystem>();
    ItemRegistry::add<Instrument>();
    ItemRegistry::add<DistantInstrument>();
    ItemRegistry::add<SEDInstrument>();
    ItemRegistry::add<FrameInstrument>();
    ItemRegistry::add<FullInstrument>();
    ItemRegistry::add<AllSkyInstrument>();
    ItemRegistry::add<HEALPixSkyInstrument>();
    ItemRegistry::add<PerspectiveInstrument>();

    // all-sky projections
    ItemRegistry::add<AllSkyProjection>();
    ItemRegistry::add<HammerAitoffProjection>();
    ItemRegistry::add<MollweideProjection>();

    // probe system and probes
    ItemRegistry::add<ProbeSystem>();
    ItemRegistry::add<Probe>();
    ItemRegistry::add<SpecialtyProbe>();
    ItemRegistry::add<SpecialtyWhenProbe>();
    ItemRegistry::add<SpecialtyWavelengthProbe>();
    ItemRegistry::add<SpecialtyWavelengthGridProbe>();
    ItemRegistry::add<SpatialGridFormProbe>();
    ItemRegistry::add<SpatialGridWhenFormProbe>();
    ItemRegistry::add<InputModelFormProbe>();
    ItemRegistry::add<ImportedSourceWeightedProbe>();
    //   .. convergence
    ItemRegistry::add<ConvergenceInfoProbe>();
    ItemRegistry::add<ConvergenceCutsProbe>();
    //   .. source
    ItemRegistry::add<LuminosityProbe>();
    ItemRegistry::add<LaunchedPacketsProbe>();
    //   .. spatial grid
    ItemRegistry::add<DensityProbe>();
    ItemRegistry::add<OpacityProbe>();
    ItemRegistry::add<TemperatureProbe>();
    ItemRegistry::add<MetallicityProbe>();
    ItemRegistry::add<VelocityProbe>();
    ItemRegistry::add<MagneticFieldProbe>();
    ItemRegistry::add<CustomStateProbe>();
    ItemRegistry::add<RadiationFieldProbe>();
    ItemRegistry::add<SecondaryDustLuminosityProbe>();
    ItemRegistry::add<SecondaryLineLuminosityProbe>();
    //   .. properties
    ItemRegistry::add<SpatialCellPropertiesProbe>();
    ItemRegistry::add<SpatialGridPlotProbe>();
    ItemRegistry::add<OpticalMaterialPropertiesProbe>();
    ItemRegistry::add<DustGrainPopulationsProbe>();
    ItemRegistry::add<DustGrainSizeDistributionProbe>();
    ItemRegistry::add<DustEmissivityProbe>();
    //   .. wavelength grid
    ItemRegistry::add<RadiationFieldWavelengthGridProbe>();
    ItemRegistry::add<DustEmissionWavelengthGridProbe>();
    ItemRegistry::add<InstrumentWavelengthGridProbe>();
    //   .. specialty
    ItemRegistry::add<DustAbsorptionPerCellProbe>();
    ItemRegistry::add<SpatialGridSourceDensityProbe>();
    ItemRegistry::add<TreeSpatialGridTopologyProbe>();
    //   .. imported source
    ItemRegistry::add<ImportedSourceLuminosityProbe>();
    ItemRegistry::add<ImportedSourceDensityProbe>();
    ItemRegistry::add<ImportedSourceMetallicityProbe>();
    ItemRegistry::add<ImportedSourceAgeProbe>();
    ItemRegistry::add<ImportedSourceVelocityProbe>();
    //   .. imported medium
    ItemRegistry::add<ImportedMediumDensityProbe>();
    ItemRegistry::add<ImportedMediumMetallicityProbe>();
    ItemRegistry::add<ImportedMediumTemperatureProbe>();
    ItemRegistry::add<ImportedMediumVelocityProbe>();

    // forms
    ItemRegistry::add<Form>();
    ItemRegistry::add<SpatialGridForm>();
    ItemRegistry::add<GenericForm>();
    ItemRegistry::add<DefaultCutsForm>();
    ItemRegistry::add<PlanarCutsForm>();
    ItemRegistry::add<PerCellForm>();
    ItemRegistry::add<LinearCutForm>();
    ItemRegistry::add<MeridionalCutForm>();
    ItemRegistry::add<AtPositionsForm>();
    ItemRegistry::add<ParallelProjectionForm>();
    ItemRegistry::add<AllSkyProjectionForm>();

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
