/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CLUMPYSPHERICALSPATIALGRID_HPP
#define CLUMPYSPHERICALSPATIALGRID_HPP

#include "SphericalClumpBVH.hpp"
#include "StructuredSphereSpatialGrid.hpp"
class Random;

////////////////////////////////////////////////////////////////////

/** An instance of the ClumpySphericalSpatialGrid class represents a spatial grid that superposes a
    classic 3D spherical grid (as offered by the StructuredSphereSpatialGrid class) with a set of
    spherical clumps. It is intended for models that need a multi-scale clump geometry -- as used,
    for example, in some models of the obscuring structure around Active Galactic Nuclei --
    combined with a background medium that need not be uniform.

    <b>%Geometry of the grid</b>

    The spatial domain of this grid is the spherical shell bounded by the inner and outer radii
    \f$r_\text{min}\ge0\f$ and \f$r_\text{max}>r_\text{min}\f$ (the inner radius may be exactly
    zero, so that the domain includes the origin). Any medium outside of this shell is ignored.

    Within this shell, a classic 3D grid in spherical coordinates is defined through three sets of
    grid points exactly as described for the StructuredSphereSpatialGrid class (radial, polar and
    azimuthal meshes). In addition, a number of spheres ("clumps") with centers and radii loaded
    from an input file are superposed on this structured grid. The clumps must be fully contained
    in the domain and are not allowed to overlap each other, but a clump \em is allowed to straddle
    the walls of one or more cells of the structured grid.

    The grid cells partitioning the spatial domain thus consist of a number of spheres (the clumps)
    plus the cells of the structured grid, each of the latter reduced in volume by the portion of
    any clumps overlapping it. As a consequence, medium properties will be uniform within each
    clump, and uniform within each structured-grid cell that is not overlapped by a clump.

    <b>Configuring the medium</b>

    During setup, as for any spatial grid, SKIRT will sample the medium properties within each
    spatial cell from the configured medium distribution. While any type of configuration will
    work, it makes most sense to specify a medium distribution that closely matches the geometry of
    this specialty grid. This can be accomplished by using a ParticleMedium with the
    UniformSmoothingKernel to match the spheres in the grid geometry. The background medium can be
    specified through any built-in geometry, such as a ShellGeometry or more likely a
    TorusGeometry. In the latter case, it is worthwile to ensure that the polar mesh configured for
    this specialty grid includes grid points exactly matching the opening angle boundaries of the
    torus.

    If the model includes multiple media types, different medium components can be combined (for
    example, one ParticleMedium for dust clumps and another one for gas clumps).

    The first four columns of the ParticleMedium import file specify the center position and radius
    of each clump. The same file can be read by this specialty grid class to define the spheres;
    the additional columns in the file are ignored. (If there are multiple ParticleMedium
    instances, the import files should be concatenated to define all spheres).

    This specialty grid class removes spheres from the imported list that (1) do not fully lie
    inside the spherical shell domain or (2) overlap any of the previously read spheres. Obviously,
    the ParticleMedium does NOT do this. As a result, the mass of these spheres will be taken into
    account while sampling the density of the spatial grid cells, most likely distorting the
    intended values. It is therefore best to ensure that the imported spheres are disjoint and
    inside the domain. */
class ClumpySphericalSpatialGrid : public StructuredSphereSpatialGrid
{
    ITEM_CONCRETE(ClumpySphericalSpatialGrid, StructuredSphereSpatialGrid,
                  "a specialty spatial grid superposing a 3D spherical grid with spherical clumps")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ClumpySphericalSpatialGrid, "Level3")

        PROPERTY_STRING(filename, "the name of the file defining the spherical clumps")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function calls the setup of the base class to build the structured grid, reads the
        input file defining the spherical clumps, and builds the data structures needed for the
        operation of the grid, including the Monte Carlo estimate of the volume of each structured
        cell that is reduced by one or more overlapping clumps. */
    void setupSelfAfter() override;

    /** The destructor destroys the BVH data structure. */
    ~ClumpySphericalSpatialGrid();

    //======================== Other Functions =======================

public:
    /** This function returns the number of cells in the grid, which equals the number of spherical
        clumps plus the number of cells in the structured grid. */
    int numCells() const override;

    /** This function returns the volume of the cell with index \f$m\f$. For a clump, this is the
        spherical volume. For a structured cell, this is the nominal cell volume minus the combined
        volume of the portions of any overlapping clumps, as estimated (and cached) during setup.
        */
    double volume(int m) const override;

    /** This function returns the diagonal of the cell with index \f$m\f$. For a clump, this is the
        diagonal of the sphere. For a structured cell, this is the distance between the cell's
        inner/lower and outer/upper corners, exactly as for the StructuredSphereSpatialGrid class.
        */
    double diagonal(int m) const override;

    /** This function returns the index \f$m\f$ of the cell that contains the position
        \f${\bf{r}}\f$, or -1 if the position is outside the domain. Cell indices \f$0\f$ through
        \f$N_\text{clumps}-1\f$ refer to the clumps, in the order in which they were retained
        during setup; cell indices \f$N_\text{clumps}\f$ and up refer to the cells of the
        structured grid, offset by \f$N_\text{clumps}\f$ from the index that
        StructuredSphereSpatialGrid::cellIndex() would assign them on its own. Because a clump can
        straddle the walls of one or more structured cells, the position is tested against the
        clumps first (using the bounding volume hierarchy, which gives an exact answer regardless
        of how thin a clump's overlap with a given structured cell is); only if it is not inside
        any clump is it located in the structured grid. */
    int cellIndex(Position bfr) const override;

    /** This function returns the central location of the cell with index \f$m\f$. For a clump,
        this is the center of the sphere. For a structured cell, the function first tries the
        nominal cell-center position; if that position happens to be inside one of the clumps, a
        random position is returned instead. */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location from the cell with index \f$m\f$. For a clump, a
        random position within the sphere is generated through analytical inversion. For a
        structured cell, a random position within the cell's nominal spherical-coordinate ranges is
        generated through analytical inversion, which is then rejected iteratively as long as it
        happens to be inside one of the clumps. */
    Position randomPositionInCell(int m) const override;

    /** This function creates and hands over ownership of a path segment generator (an instance of
        a PathSegmentGenerator subclass) appropriate for this spatial grid type. */
    std::unique_ptr<PathSegmentGenerator> createPathSegmentGenerator() const override;

protected:
    /** This function writes the intersection of the grid structure with the xy plane to the
        specified SpatialGridPlotFile object. */
    void write_xy(SpatialGridPlotFile* outfile) const override;

    /** This function writes the intersection of the grid structure with the xz plane to the
        specified SpatialGridPlotFile object. */
    void write_xz(SpatialGridPlotFile* outfile) const override;

    /** This function writes the intersection of the grid structure with the yz plane to the
        specified SpatialGridPlotFile object. */
    void write_yz(SpatialGridPlotFile* outfile) const override;

    /** This function writes 3D information for the grid structure to the specified
        SpatialGridPlotFile object. It calls the base class implementation to write the structured
        grid, and then adds a coarse wireframe sphere for each clump. */
    void write_xyz(SpatialGridPlotFile* outfile) const override;

    /** This function returns a more conservative (larger) epsilon scale factor than the base
        class's default, since paths in this class cross more boundaries of more varied kinds
        (structured-cell walls interleaved with clump entries and exits) within a given region. */
    double meshEpsilonScale() const override;

    //======================== Data Members ========================

private:
    // convenience alias for the clump type defined and used by SphericalClumpBVH
    using Clump = SphericalClumpBVH::Clump;

    // array defining the clumps
    int _numClumps{0};      // the number of clumps AND the offset of the first structured cell index
    vector<Clump> _clumps;  // index on m, assuming m < _numClumps

    // the bounding volume hierarchy that allows efficient querying over the clumps
    SphericalClumpBVH* _bvh{nullptr};

    // allow our path segment generator to access our private data members
    class MySegmentGenerator;
    friend class MySegmentGenerator;
};

////////////////////////////////////////////////////////////////////

#endif
