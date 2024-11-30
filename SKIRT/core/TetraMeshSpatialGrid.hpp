/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TETRAMESHSPATIALGRID_HPP
#define TETRAMESHSPATIALGRID_HPP

#include "BoxSpatialGrid.hpp"
#include "DensityInCellInterface.hpp"
#include "Log.hpp"
#include "Medium.hpp"
#include "PathSegmentGenerator.hpp"

class tetgenio;

//////////////////////////////////////////////////////////////////////

/** TetraMeshSpatialGrid is a concrete subclass of the SpatialGrid class. It represents a
    three-dimensional grid based on a Tetra tesselation of the cuboidal spatial domain of the
    simulation. See the TetraMeshSnapshot class for more information on Tetra tesselations.

    The class offers several options for determining the positions of the sites generating the
    Tetra tesselation. A specified number of sites can be distributed randomly over the domain,
    either uniformly or with the same overall density distribution as the medium. Alternatively,
    the positions can be copied from the sites in the imported distribution(s).

    Furthermore, the user can opt to perform a relaxation step on the site positions to avoid
    overly elongated cells. */
class TetraMeshSpatialGrid : public BoxSpatialGrid
{
    /** The enumeration type indicating the policy for determining the positions of the sites. */
    ENUM_DEF(Policy, Uniform, CentralPeak, DustDensity, ElectronDensity, GasDensity, File, ImportedSites, ImportedMesh)
        ENUM_VAL(Policy, Uniform, "random from uniform distribution")
        ENUM_VAL(Policy, CentralPeak, "random from distribution with a steep central peak")
        ENUM_VAL(Policy, DustDensity, "random from dust density distribution")
        ENUM_VAL(Policy, ElectronDensity, "random from electron density distribution")
        ENUM_VAL(Policy, GasDensity, "random from gas density distribution")
        ENUM_VAL(Policy, ImportedSites, "positions of particles, sites or cells in imported distribution")
    ENUM_END()

    ITEM_CONCRETE(TetraMeshSpatialGrid, BoxSpatialGrid, "a Tetra tessellation-based spatial grid")
        ATTRIBUTE_TYPE_DISPLAYED_IF(TetraMeshSpatialGrid, "Level2")

        PROPERTY_ENUM(policy, Policy, "the policy for determining the positions of the sites")
        ATTRIBUTE_DEFAULT_VALUE(policy, "DustDensity")

        PROPERTY_INT(numSites, "the number of random sites to be used as vertices")
        ATTRIBUTE_MIN_VALUE(numSites, "4")
        ATTRIBUTE_DEFAULT_VALUE(numSites, "500")
        ATTRIBUTE_RELEVANT_IF(numSites, "policyUniform|policyCentralPeak|policyDustDensity|"
                                        "policyElectronDensity|policyGasDensity")

        PROPERTY_BOOL(refine, "refine the grid to have higher quality cells by adding more vertices")
        ATTRIBUTE_DEFAULT_VALUE(refine, "false");

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** The destructor releases the Tetra mesh if this object owns it. */
    ~TetraMeshSpatialGrid();

protected:
    /** This function verifies that the attributes have been appropriately set, generates or
        retrieves the site positions for constructing the Tetra tessellation according to the
        configured policy, and finally constructs the Tetra tessellation through an instance of
        the TetraMeshSnapshot class. */
    void setupSelfBefore() override;

    //==================== Private construction ====================
private:
    class Node;
    class Face;
    class Tetra;

    void buildMesh();

    void buildDelaunay(tetgenio& out);

    void refineDelaunay(tetgenio& in, tetgenio& out);

    void storeTetrahedra(const tetgenio& out, bool storeVertices);

    /** Private function to recursively build a binary search tree (see
        en.wikipedia.org/wiki/Kd-tree) */
    Node* buildTree(vector<int>::iterator first, vector<int>::iterator last, int depth) const;

    /** This private function builds data structures that allow accelerating the operation of the
        cellIndex() function. It assumes that the Voronoi mesh has already been built.

        The domain is partitioned using a linear cubodial grid into cells that are called \em
        blocks. For each block, the function builds and stores a list of all Voronoi cells that
        possibly overlap that block. Retrieving the list of cells possibly overlapping a given
        point in the domain then comes down to locating the block containing the point (which is
        trivial since the grid is linear). The current implementation uses a Voronoi cell's
        enclosing cuboid to test for intersection with a block. Performing a precise intersection
        test is \em really slow and the shortened block lists don't substantially accelerate the
        cellIndex() function.

        To further reduce the search time within blocks that overlap with a large number of cells,
        the function builds a binary search tree on the cell sites for those blocks (see for example
        <a href="http://en.wikipedia.org/wiki/Kd-tree">en.wikipedia.org/wiki/Kd-tree</a>). */
    void buildSearchPerBlock();

    //======================== Other Functions =======================

public:
    /** This function returns the number of cells in the grid. */
    int numCells() const override;

    /** This function returns the volume of the cell with index \f$m\f$. */
    double volume(int m) const override;

    /** This function returns the approximate diagonal of the cell with index \f$m\f$. For a
        Tetra grid, it returns \f$(3V)^(1/3)\f$ where \f$V\f$ is the volume of the cell. This
        corresponds to the correct diagonal only for cubical cells. */
    double diagonal(int m) const override;

    /** This function returns the index of the cell that contains the position \f${\bf{r}}\f$. */
    int cellIndex(Position bfr) const override;

    /** This function returns the central location of the cell with index \f$m\f$. In this class
        the function returns the centroid of the Tetra cell. */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location from the cell with index \f$m\f$. */
    Position randomPositionInCell(int m) const override;

    /** This function creates and hands over ownership of a path segment generator (an instance of
        a PathSegmentGenerator subclass) appropriate for this spatial grid type. For the Tetra
        mesh grid, the path segment generator is actually implemented in the TetraMeshSnapshot
        class. */
    std::unique_ptr<PathSegmentGenerator> createPathSegmentGenerator() const override;

    /** This function outputs the grid plot files; it is provided here because the regular
        mechanism does not apply. The function reconstructs the Tetra tesselation in order to
        produce the coordinates of the Tetra cell vertices. */
    void writeGridPlotFiles(const SimulationItem* probe) const override;

    //======================== Data Members ========================

private:
    // data members initialized by setupSelfBefore()
    Log* _log{nullptr};

    // data members initialized during configuration
    double _eps{0.};  // small fraction of extent

    int _numTetra;
    int _numVertices;  // vertices are added/removed as the grid is built and refined
    vector<Tetra*> _tetrahedra;
    vector<Vec*> _vertices;
    vector<Vec*> _centroids;

    // data members initialized by BuildSearch()
    int _nb{0};                       // number of blocks in each dimension (limit for indices i,j,k)
    int _nb2{0};                      // nb*nb
    int _nb3{0};                      // nb*nb*nb
    vector<vector<int>> _blocklists;  // list of cell indices per block, indexed on i*_nb2+j*_nb+k
    vector<Node*> _blocktrees;        // root node of search tree or null for each block, indexed on i*_nb2+j*_nb+k

    // allow our path segment generator to access our private data members
    class MySegmentGenerator;
    friend class MySegmentGenerator;
};

//////////////////////////////////////////////////////////////////////

#endif
