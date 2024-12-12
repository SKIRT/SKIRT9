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
#include <array>

class tetgenio;

//////////////////////////////////////////////////////////////////////

/** TetraMeshSpatialGrid is a concrete subclass of the SpatialGrid class. It represents a
    three-dimensional grid based on a tetrahedralization of the cuboidal spatial domain of the
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
    ENUM_DEF(Policy, Uniform, CentralPeak, DustDensity, ElectronDensity, GasDensity, File)
        ENUM_VAL(Policy, Uniform, "random from uniform distribution")
        ENUM_VAL(Policy, CentralPeak, "random from distribution with a steep central peak")
        ENUM_VAL(Policy, DustDensity, "random from dust density distribution")
        ENUM_VAL(Policy, ElectronDensity, "random from electron density distribution")
        ENUM_VAL(Policy, GasDensity, "random from gas density distribution")
        ENUM_VAL(Policy, File, "loaded from text column data file")
    ENUM_END()

    ITEM_CONCRETE(TetraMeshSpatialGrid, BoxSpatialGrid, "a tetrahedral spatial grid")
        ATTRIBUTE_TYPE_DISPLAYED_IF(TetraMeshSpatialGrid, "Level2")

        PROPERTY_ENUM(policy, Policy, "the policy for determining the positions of the sites")
        ATTRIBUTE_DEFAULT_VALUE(policy, "DustDensity")

        PROPERTY_INT(numSites, "the number of random sites to be used as vertices")
        ATTRIBUTE_MIN_VALUE(numSites, "4")
        ATTRIBUTE_DEFAULT_VALUE(numSites, "500")
        ATTRIBUTE_RELEVANT_IF(numSites, "policyUniform|policyCentralPeak|policyDustDensity|"
                                        "policyElectronDensity|policyGasDensity")

        PROPERTY_STRING(filename, "the name of the file containing the site positions")
        ATTRIBUTE_RELEVANT_IF(filename, "policyFile")

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
    /** Private class to hold all infromation for a face of a tetrahedron to allow photon
        traversal through it. */
    struct Face
    {
        int _ntetra;  // index of neighbouring tetrahedron
        int _nface;   // neighbouring face index
        Vec _normal;  // outward facing normal

        Face(){};

        Face(int ntetra, int nface, Vec normal) : _ntetra(ntetra), _nface(nface), _normal(normal) {}
    };

    /** Private class to hold all information for a tetrahedron to fully describe its geometry
        and allow photon traversal through it. */
    class Tetra
    {
    private:
        const vector<Vec>& _vertices;  // reference to all vertices, could be passed as arg but this is more convenient
        Box _extent;                   // bounding box of the tetrahedron
        Vec _centroid;                 // barycenter of the tetrahedron
        std::array<int, 4> _vertexIndices;  // indices of the vertices in the full list
        std::array<Face, 4> _faces;         // face information

    public:
        Tetra(const vector<Vec>& vertices, const std::array<int, 4>& vertexIndices, const std::array<Face, 4>& faces);

        int findEnteringFace(const Vec& pos, const Direction& dir) const;

        bool contains(const Position& bfr) const;

        double generateBarycentric(double& s, double& t, double& u) const;

        Position generatePosition(Random* random) const;

        Vec vertex(int t) const;

        Vec edge(int t1, int t2) const;

        double volume() const;

        double diagonal() const;

        const std::array<Face, 4>& faces() const;

        const Vec& centroid() const;

        const Box& extent() const;
    };

    class BlockGrid;

    /** This private function removes vertices that are outside the domain. */
    void removeOutside();

    /** This private function builds the tetrahedral mesh. Starts by constructing the Delaunay 
        tetrahedralization and optionally refines it if the refine option is set to true. */
    void buildMesh();

    /** This private function constructs the Delaunay tetrahedralization using the _vertices. */
    void buildDelaunay(tetgenio& out);

    /** This private function refines the Delaunay tetrahedralization to improve cell quality.
        The refinement process is controlled by TetGen with default quality parameters. The
        input is the initial Delaunay tetrahedralization in a tetgenio container. The output
        is a tetgenio container with the refined Delaunay tetrahedralization. */
    void refineDelaunay(tetgenio& in, tetgenio& out);

    /** This private function stores the tetrahedra and vertices from the final tetgenio 
        container into the class members. The input is the tetgenio container with the
        final tetrahedralization. The storeVertices parameter indicates whether to
        overwrite the vertices with those from the tetgenio container. This function
        also logs some cell statistics after storing the data. */
    void storeTetrahedra(const tetgenio& final, bool storeVertices);

    /** This private function builds the search data structure for the tetrahedral mesh.
        See the private BlockGrid class for more information. */
    void buildSearch();

    //======================= Interrogation =======================
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

    //===================== Path construction =====================

    /** This function creates and hands over ownership of a path segment generator (an instance of
        a PathSegmentGenerator subclass) appropriate for this spatial grid type. For the Tetra
        mesh grid, the path segment generator is actually implemented in the TetraMeshSnapshot
        class. */
    std::unique_ptr<PathSegmentGenerator> createPathSegmentGenerator() const override;

    //===================== Output =====================

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

    int _numCells;     // number of Tetra cells and centroids
    int _numVertices;  // vertices are added/removed as the grid is built and refined
    vector<Tetra> _tetrahedra;
    vector<Vec> _vertices;

    // smart grid that contains all cell information
    // allows for efficiently locating the cell at a given location
    BlockGrid* _blocks{nullptr};

    // allow our path segment generator to access our private data members
    class MySegmentGenerator;
    friend class MySegmentGenerator;
};

//////////////////////////////////////////////////////////////////////

#endif
