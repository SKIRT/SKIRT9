/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CLUMPYSPHERICALSPATIALGRID_HPP
#define CLUMPYSPHERICALSPATIALGRID_HPP

#include "Mesh.hpp"
#include "SphereSpatialGrid.hpp"
class Random;

////////////////////////////////////////////////////////////////////

/** An instance of the ClumpySphericalSpatialGrid class represents a spatial grid that superposes a
    classic 3D spherical grid (as offered by the Sphere3DSpatialGrid class) with a set of spherical
    clumps (as offered by the ClumpyTorusSpatialGrid class). It is intended for models that need
    the multi-scale clump geometry used for AGN torus models, but without the torus's conical
    opening-angle restriction, and/or with a background medium that is not itself uniform (unlike
    the single background cell of ClumpyTorusSpatialGrid).

    <b>%Geometry of the grid</b>

    The spatial domain of this grid is the spherical shell bounded by the inner and outer radii
    \f$r_\text{min}\ge0\f$ and \f$r_\text{max}>r_\text{min}\f$ (unlike ClumpyTorusSpatialGrid, the
    inner radius may be exactly zero, so that the domain includes the origin). Any medium outside
    of this shell is ignored.

    Within this shell, a classic 3D grid in spherical coordinates is defined through three sets of
    grid points exactly as described for the Sphere3DSpatialGrid class (radial, polar and azimuthal
    meshes). In addition, a number of spheres ("clumps") with centers and radii loaded from an input
    file are superposed on this structured grid. The clumps must be fully contained in the domain
    and are not allowed to overlap each other, but -- in contrast to ClumpyTorusSpatialGrid -- a
    clump \em is allowed to straddle the walls of one or more cells of the structured grid.

    The grid cells partitioning the spatial domain thus consist of a number of spheres (the clumps)
    plus the cells of the structured grid, each of the latter reduced in volume by the portion of
    any clumps overlapping it. As a consequence, medium properties will be uniform within each
    clump, and (given a uniform density field) uniform within each structured-grid cell that is not
    overlapped by a clump.

    <b>Configuring the medium</b>

    As for the ClumpyTorusSpatialGrid class, the first four columns of the ParticleMedium import
    file specify the center position and radius of each clump, and the same file can be read by
    this specialty grid class to define the spheres; see the ClumpyTorusSpatialGrid documentation
    for further details and caveats regarding matching the medium configuration to the grid
    geometry. */
class ClumpySphericalSpatialGrid : public SphereSpatialGrid
{
    ITEM_CONCRETE(ClumpySphericalSpatialGrid, SphereSpatialGrid,
                  "a specialty spatial grid superposing a 3D spherical grid with spherical clumps")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ClumpySphericalSpatialGrid, "Level3")

        PROPERTY_ITEM(meshRadial, Mesh, "the bin distribution in the radial direction")
        ATTRIBUTE_DEFAULT_VALUE(meshRadial, "LinMesh")

        PROPERTY_ITEM(meshPolar, Mesh, "the bin distribution in the polar direction")
        ATTRIBUTE_DEFAULT_VALUE(meshPolar, "LinMesh")

        PROPERTY_ITEM(meshAzimuthal, Mesh, "the bin distribution in the azimuthal direction")
        ATTRIBUTE_DEFAULT_VALUE(meshAzimuthal, "LinMesh")

        PROPERTY_STRING(filename, "the name of the file defining the spherical clumps")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function sets up the structured grid (as for Sphere3DSpatialGrid), reads the input file
        defining the spherical clumps, and builds the data structures needed for the operation of
        the grid, including the Monte Carlo estimate of the volume of each structured cell that is
        reduced by one or more overlapping clumps. */
    void setupSelfAfter() override;

    /** The destructor destroys the custom BVH data structure. */
    ~ClumpySphericalSpatialGrid();

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the grid, which is 3 for this class. */
    int dimension() const override;

    /** This function returns the number of cells in the grid, which equals the number of spherical
        clumps plus the number of cells in the structured grid. */
    int numCells() const override;

    /** This function returns the volume of the cell with index \f$m\f$. For a clump, this is the
        spherical volume. For a structured cell, this is the nominal cell volume minus the combined
        volume of the portions of any overlapping clumps, as estimated (and cached) during setup. */
    double volume(int m) const override;

    /** This function returns the diagonal of the cell with index \f$m\f$. For a clump, this is the
        diagonal of the sphere. For a structured cell, this is the distance between the cell's
        inner/lower and outer/upper corners, exactly as for the Sphere3DSpatialGrid class. */
    double diagonal(int m) const override;

    /** This function returns the index \f$m\f$ of the cell that contains the position
        \f${\bf{r}}\f$. */
    int cellIndex(Position bfr) const override;

    /** This function returns the central location of the cell with index \f$m\f$. For a clump,
        this is the center of the sphere. For a structured cell, the function first tries the
        nominal cell-center position; if that position happens to be inside one of the clumps (as
        determined through an exact query of the bounding volume hierarchy), a random position is
        returned instead. */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location from the cell with index \f$m\f$. For a clump, a
        random position within the sphere is generated through analytical inversion. For a
        structured cell, a random position within the cell's nominal spherical-coordinate ranges is
        generated through analytical inversion, which is then rejected iteratively -- based on an
        exact query of the bounding volume hierarchy -- as long as it happens to be inside one of
        the clumps. */
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

    /** This function is intentionally left unimplemented for this grid type. */
    void write_xyz(SpatialGridPlotFile* outfile) const override;

    //======================== Helper Functions =======================

private:
    /** This private function returns the local structured-cell index (i.e. not counting the clump
        cells) corresponding to the radial index \f$i\f$, the polar index \f$j\f$, and azimuthal
        index \f$k\f$, using the same correspondence as for Sphere3DSpatialGrid. */
    int structuredIndex(int i, int j, int k) const;

    /** This function obtains the spherical coordinates for the corners of the structured cell with
        local index \f$s\f$ (i.e. not counting the clump cells), exactly as for Sphere3DSpatialGrid.
        If the index is out of range, the function returns false and leaves the provided arguments
        unchanged. */
    bool getCoords(int s, double& rmin, double& thetamin, double& phimin, double& rmax, double& thetamax,
                   double& phimax) const;

    //======================== Data Members ========================

private:
    // data type defining a single clump
    class Clump
    {
        double _x, _y, _z, _r;
        int _numOverlappingCells{0};

    public:
        Clump(double x, double y, double z, double r) : _x{x}, _y{y}, _z{z}, _r{r} {}
        Position center() const { return Position(_x, _y, _z); }
        double radius() const { return _r; }
        Box bounds() const { return Box(_x - _r, _y - _r, _z - _r, _x + _r, _y + _r, _z + _r); }
        void setNumOverlappingCells(int n) { _numOverlappingCells = n; }
        int numOverlappingCells() const { return _numOverlappingCells; }
    };

    // array defining the clumps
    int _numClumps{0};      // the number of clumps AND the offset of the first structured cell index
    vector<Clump> _clumps;  // index on m, assuming m < _numClumps

    // the custom bounding volume hierarchy that allows efficient querying for our purposes
    class BVH;
    BVH* _bvh{nullptr};

    // the structured grid, exactly as for Sphere3DSpatialGrid
    int _Nr{0};
    int _Ntheta{0};
    int _Nphi{0};
    int _Ncells{0};
    Array _rv;
    Array _thetav;
    Array _phiv;
    Array _cv;
    Array _sinv;
    Array _cosv;

    // the volume of each structured cell (local index s, 0 <= s < _Ncells), precomputed and cached
    // during setup as the nominal cell volume minus the volume of any overlapping clump portions
    vector<double> _cellVolume;

    // small value used to progress a path
    double _eps{0.};

    // allow our path segment generator to access our private data members
    class MySegmentGenerator;
    friend class MySegmentGenerator;
};

////////////////////////////////////////////////////////////////////

#endif
