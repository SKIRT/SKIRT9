/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TETRAMESHSPATIALGRID_HPP
#define TETRAMESHSPATIALGRID_HPP

#include "BoxSpatialGrid.hpp"
#include "Log.hpp"
#include "PathSegmentGenerator.hpp"
#include <array>

class tetgenio;

//////////////////////////////////////////////////////////////////////

/** TetraMeshSpatialGrid is a concrete subclass of SpatialGrid representing a 3D tetrahedral mesh
    generated from a set of vertices. The resulting grid fully covers the domain by including the
    8 corners of the domain as vertices. Since a tetrahedral mesh always fills the convex hull of
    the vertices, the domain is always fully covered. This grid will thus always have at least 8
    vertices. The grid is constructed using the open-source library TetGen version 1.6.0 (released
    on August 31, 2020). TetGen is an advanced C++ tetrahedral mesh generator with many features.
    This class uses the following key features of TetGen:

    - Delaunay Tetrahedralization:
      Generate a unique Delaunay tetrahedral mesh from a given set of vertices.

    - %Mesh Refinement:
      Refine the tetrahedral mesh using TetGen's Delaunay refinement algorithm. This option
      is enabled through the \em refine property of this class.

    It should be noted that TetGen is a single-threaded library, but the algorithms used are
    generally quite fast. The refinement process is by far the most time-consuming part of the
    mesh generation.

    The 3D Delaunay tetrahedralization often contains cells less suited for a computational grid.
    While the 2D Delaunay triangulation maximizes the smallest angle in each triangle, resulting
    in more regular cells, the 3D version lacks this property. Consequently, 3D tetrahedralizations
    frequently contain elongated, irregular cells. To address this, Delaunay refinement algorithms
    can be used. TetGen implements its own refinement process, applying local mesh operations like
    vertex smoothing, edge and face swapping, edge contraction, and vertex insertion. This
    refinement aims to optimize quality metrics such as the radius-edge ratio and the minimum
    dihedral angle. All refinement parameters are kept at their default values provided by TetGen.
    However, it is uncertain whether this refinement process will greatly improve the grids for
    radiative transfer applications.

    The positions of the vertices used for generating the tetrahedral mesh are determined by the
    \em policy property. The available policies are:

    - Uniform: Randomly sampled from a uniform distribution.
    - CentralPeak: Randomly sampled from a distribution with a steep central peak.
    - DustDensity: Randomly sampled based on the dust density distribution.
    - ElectronDensity: Randomly sampled based on the electron density distribution.
    - GasDensity: Randomly sampled based on the gas density distribution.
    - File: Loaded from a column data file specified by the \em filename property,
            containing vertex coordinates (x, y, z) in each column.

    Vertices are removed if they lie outside the simulation domain or are too close to another.
*/
class TetraMeshSpatialGrid : public BoxSpatialGrid
{
    /** The enumeration type indicating the policy for determining the positions of the vertices. */
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

        PROPERTY_ENUM(policy, Policy, "the policy for determining the positions of the vertices")
        ATTRIBUTE_DEFAULT_VALUE(policy, "DustDensity")

        PROPERTY_INT(numSamples, "the number of random positions to be used as vertices")
        ATTRIBUTE_MIN_VALUE(numSamples, "4")
        ATTRIBUTE_DEFAULT_VALUE(numSamples, "500")
        ATTRIBUTE_RELEVANT_IF(numSamples, "policyUniform|policyCentralPeak|policyDustDensity|"
                                          "policyElectronDensity|policyGasDensity")

        PROPERTY_STRING(filename, "the name of the file containing the vertex positions")
        ATTRIBUTE_RELEVANT_IF(filename, "policyFile")

        PROPERTY_BOOL(refine, "refine the grid by performing local mesh operations")
        ATTRIBUTE_DEFAULT_VALUE(refine, "false");

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This destructor releases the BlockGrid search structure. */
    ~TetraMeshSpatialGrid();

protected:
    /** This function verifies that the attributes are correctly set, generates or retrieves
        vertex positions based on the configured policy, builds (and optionally refines) the
        tetrahedralization, and constructs the search structure to optimize the \em CellIndex
        function. */
    void setupSelfBefore() override;

    //==================== Private data types ====================

private:
    /** Private struct that represents a face of a tetrahedron, storing all relevant
        information for photon traversal. */
    struct Face
    {
        int _ntetra;  // cell index of neighbouring tetrahedron
        int _nface;   // index of the equivalent face in the neighbouring tetrahedron [0, 3]
        Vec _normal;  // outward facing normal

        Face() : _ntetra(0), _nface(0) {}

        Face(int ntetra, int nface, Vec normal) : _ntetra(ntetra), _nface(nface), _normal(normal) {}
    };

    /** Alias for a fixed array of 4 integer indices. */
    using FourIndices = std::array<int, 4>;

    /** Alias for a fixed array of 4 Faces. */
    using FourFaces = std::array<Face, 4>;

    /** Private class that represents a tetrahedron, storing all relevant information to fully
        describe its geometry and allow photon traversal through it. */
    class Tetra
    {
    private:
        const vector<Vec>& _vertices;  // reference to all vertices, could be passed as arg but this is more convenient
        Box _extent;                   // bounding box of the tetrahedron
        Vec _centroid;                 // barycenter of the tetrahedron
        FourIndices _vertexIndices;    // indices of the vertices in the global vertex list
        FourFaces _faces;              // faces of the tetrahedron, see struct Face

    public:
        /** This constructor initializes the tetrahedron by setting its fields and calculating the
            bounding box and centroid. */
        Tetra(const vector<Vec>& vertices, const FourIndices& vertexIndices, const FourFaces& faces);

        /** This function finds a face that is not the leaving face which can act as the entering
            face in the traversal algorithm. */
        int findEnteringFace(const Vec& pos, const Direction& dir) const;

        /** This function checks if the given position is contained inside the tetrahedron.
            It first checks if the position is inside the bounding box of the tetrahedron. */
        bool contains(const Position& bfr) const;

        /** This function generates three random barycentric coordinates for uniformly sampling
            inside this tetrahedron. The fourth coordinate is calculated by ensuring their sum
            equals 1, i.e. r=1-s-t-u.
            Source: Generating Random Points in a Tetrahedron: DOI 10.1080/10867651.2000.10487528 */
        double generateBarycentric(double& s, double& t, double& u) const;

        /** This function generates a random position inside the tetrahedron by generating random
             barycentric coordinates and using the vertices of the tetrahedron. */
        Position generatePosition(Random* random) const;

        /** This function returns the vertex with index \f$t\f$ of the tetrahedron, where
            \f$t \in \{0, 1, 2, 3\}\f$. */
        Vec vertex(int t) const;

        /** This function returns the edge from vertex \f$t1\f$ to \f$t2\f$ of the tetrahedron,
            where \f$t1, t2 \in \{0, 1, 2, 3\}\f$. */
        Vec edge(int t1, int t2) const;

        /** This function calculates and returns the volume of the tetrahedron. */
        double volume() const;

        /** This function calculates and returns the approximate diagonal of the tetrahedron.
            It calculates the square root of the average of the squared edge lengths. */
        double diagonal() const;

        /** This function returns a reference to an array of the faces of the tetrahedron. */
        const FourFaces& faces() const;

        /** This function returns the centroid of the tetrahedron. */
        const Vec& centroid() const;

        /** This function returns the extent of the tetrahedron. */
        const Box& extent() const;
    };

    //==================== Private construction ====================

private:
    /** This private helper class organizes the cells into cuboidal blocks in a smart grid, such
        that it is easy to retrieve all cells inside a certain block given a position. */
    class BlockGrid;

    /** This private function removes vertices that are outside the domain or too close to other vertices. */
    void removeInvalid();

    /** This private function adds the 8 corners of the domain to the vertex list. This way the full
        domain will be tetrahedralized. */
    void addCorners();

    /** This private function builds the tetrahedral mesh. It starts by constructing the Delaunay
        tetrahedralization and optionally refines it if the \em refine option is set to true. */
    void buildMesh();

    /** This private function constructs the Delaunay tetrahedralization using the vertices obtained
        from the \em policy. The output is placed inside the tetgenio reference. */
    void buildDelaunay(tetgenio& out);

    /** This private function refines the Delaunay tetrahedralization to improve cell quality. The
        refinement process is controlled by TetGen with default quality parameters. The input is
        the initial Delaunay tetrahedralization in the \em in tetgenio reference. The output, with
        the refined Delaunay tetrahedralization, is placed inside the \em out tetgenio reference. */
    void refineDelaunay(tetgenio& in, tetgenio& out);

    /** This private function stores the tetrahedra and vertices from the \em final tetgenio container
        into the \em TetraMeshSpatialGrid members. The input is the tetgenio reference with the final
        tetrahedralization. The \em storeVertices parameter indicates whether to overwrite the vertices
        with those from the tetgenio container. This function also logs some cell statistics after
        it finishes transferring the data. */
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
        tetrahedron, it takes the square root of the average of the squared edge lengths. */
    double diagonal(int m) const override;

    /** This function returns the index of the cell that contains the position \f${\bf{r}}\f$.
        The function uses the data structure stored in the \em BlockGrid to accelerate the
        search. If no cell is found to contain this position, the function returns -1. */
    int cellIndex(Position bfr) const override;

    /** This function returns the centroid of the tetrahedron with index \f$m\f$. */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location from the tetrahedron with index \f$m\f$. */
    Position randomPositionInCell(int m) const override;

    //===================== Path construction =====================

public:
    /** This function creates and returns ownership of a path segment generator suitable for the
        tetrahedral spatial grid, implemented as a private \em PathSegmentGenerator subclass. The
        algorithm for constructing the path is taken from Maria et al. (2017).

        The algorithm uses Plücker coordinates to identify the exit face of a ray inside a given
        tetrahedron. Plücker coordinates are a set of six values that describe a directed line in 3D space:
        \f[\mathbf{\pi}_R = \{\mathbf{k} : \mathbf{k} \times \mathbf{r}\} = \{\mathbf{U}_R : \mathbf{V}_R\}\f]
        where \f$\mathbf{r}\f$ is a position along the ray, and \f$\mathbf{k}\f$ is the ray's direction.
        The Plücker product is defined as:
        \f[\mathbf{\pi}_R \odot \mathbf{\pi}_S = \mathbf{U}_R \cdot \mathbf{V}_S + \mathbf{U}_S \cdot \mathbf{V}_R\f]
        This product determines the relative orientation between two rays:
        \f[\mathbf{\pi}_R \odot \mathbf{\pi}_S = \begin{cases}
        > 0 & \iff S \text{ goes counterclockwise around } R\\
        < 0 & \iff S \text{ goes clockwise around } R\\
        = 0 & \iff S \text{ intersects or is parallel to } R \end{cases} \f]
        A ray exits a tetrahedron through a particular face if the Plücker products with all three
        clockwise-ordered edges of that face are negative. The algorithm in Maria et al. (2017) optimizes
        this by requiring only two Plücker products to be evaluated. Our implementation is described below.

        In the first step, the function checks whether the start point is inside the domain. If so, the current
        point is simply initialized to the start point. If not, the function computes the path segment to the
        first intersection with one of the domain walls and moves the current point inside the domain. Next,
        the function determines the cell index of the tetrahedron containing the current point. If none is
        found, the path is terminated. Before the traversal algorithm can commence, a non-leaving face must
        be identified. This face acts as the entry face for the ray. Note that this face does not necessarily
        have to be the actual entry face. This task is handled by the \em findEnteringFace function of the
        \em Tetra class.

        Next, the traversal algorithm begins. The entering face is labeled as face 0, with its opposing vertex
        labeled as vertex 0. We start by evaluating the Plücker product of the ray with the edge \f$1 \to 0\f$.
        Based on whether this product is positive or negative, the next Plücker product to evaluate is determined.
        If the product is negative (i.e., clockwise orientation), we evaluate the product with the edge in the
        clockwise direction viewed from vertex 0. The same applies for a positive product. The exit face is
        determined using the decision tree described in Maria et al. (2017). The distance traveled through the
        cell is calculated using a simple line-plane intersection:
        \f[s_i = \frac{\mathbf{n} \cdot (\mathbf{v} - \mathbf{r})}{\mathbf{n} \cdot \mathbf{k}}\f]
        where \f$\mathbf{n}\f$ is the outward-pointing normal of the face, \f$\mathbf{v}\f$ is any vertex on
        the exit face, and \f$\mathbf{r}\f$ is the current position.

        Although the algorithm described in Maria et al. (2017) works even if one of the Plücker products is zero,
        we revert to a plane intersection algorithm in such cases. This approach is similar to the one used in the
        \em VoronoiMeshSnapshot class, where the closest intersection distance with all faces is found.

        The algorithm continues until the exit face lies on the convex hull boundary. At this point, the path is
        terminated. If the exit face is not found, which should only rarely happen due to computational
        inaccuracies, the current point is advanced by a small distance, and the cell index is recalculated. */
    std::unique_ptr<PathSegmentGenerator> createPathSegmentGenerator() const override;

    //===================== Output =====================

public:
    /** This function outputs the grid plot files. It writes each tetrahedral face as a triangle
        to the grid plot file. */
    void writeGridPlotFiles(const SimulationItem* probe) const override;

    //======================== Data Members ========================

private:
    // data members initialized by setupSelfBefore()
    Log* _log{nullptr};
    double _eps{0.};  // small fraction of extent

    // data members describing the tetrahedralization
    int _numCells;     // total number of tetrahedra
    int _numVertices;  // total number of vertices
    vector<Tetra> _tetrahedra;
    vector<Vec> _vertices;

    // smart grid that organizes the tetrahedra into blocks
    BlockGrid* _blocks{nullptr};

    // allow our path segment generator to access our private data members
    class MySegmentGenerator;
    friend class MySegmentGenerator;
};

//////////////////////////////////////////////////////////////////////

#endif
