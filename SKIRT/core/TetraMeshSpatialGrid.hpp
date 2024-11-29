/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TETRAMESHSPATIALGRID_HPP
#define TETRAMESHSPATIALGRID_HPP

#include "BoxSpatialGrid.hpp"
#include "DensityInCellInterface.hpp"
#include "Medium.hpp"

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
class TetraMeshSpatialGrid : public BoxSpatialGrid, public DensityInCellInterface
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
        ATTRIBUTE_MIN_VALUE(numSites, "5")
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

private:
    class Node;
    class Face;
    class Tetra;

    /** Build the Delaunay triangulation */
    void buildDelaunay(tetgenio& out);

    void refineDelaunay(tetgenio& in, tetgenio& out);

    void buildMesh();

    void storeTetrahedra(const tetgenio& out, bool restoreVertices);

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

    /** This function implements the DensityInCellInterface interface. It returns the number
        density for medium component \f$h\f$ in the grid cell with index \f$m\f$. For a Tetra
        mesh grid, this interface can be offered only if the "ImportedMesh" policy has been
        configured, and the medium system consists of a single component, which then, by
        definition, supplies the Tetra mesh. Thus, the component index \f$h\f$ passed to this
        function should always be zero; in fact, its value is actually ignored. */
    double numberDensity(int h, int m) const override;

protected:
    /** This function is used by the interface() function to ensure that the receiving item can
        actually offer the specified interface. If the requested interface is the
        DensityInCellInterface, the implementation in this class returns true if the "ImportedMesh"
        policy has been configured and the medium system consists of a single component, and false
        otherwise. For other requested interfaces, the function invokes its counterpart in the base
        class. */
    bool offersInterface(const std::type_info& interfaceTypeInfo) const override;

private:
    Node* buildTree(vector<int>::iterator first, vector<int>::iterator last, int depth) const;

    void buildSearchPerBlock();

    //======================== Data Members ========================

private:
    Log* _log{nullptr};

    // data members initialized during configuration
    double _eps{0.};  // small fraction of extent
    int _numTetra;
    int _numVertices; // vertices are added/removed as the grid is built and refined

    // data members initialized when processing snapshot input and further completed by BuildMesh()
    vector<Tetra*> _tetrahedra;
    vector<Vec*> _vertices;
    vector<Vec*> _centroids;

    // data members initialized when processing snapshot input, but only if a density policy has been set
    Array _rhov;       // density for each cell (not normalized)
    double _mass{0.};  // total effective mass

    // data members initialized by BuildSearch()
    int _nb{0};                       // number of blocks in each dimension (limit for indices i,j,k)
    int _nb2{0};                      // nb*nb
    int _nb3{0};                      // nb*nb*nb
    vector<vector<int>> _blocklists;  // list of cell indices per block, indexed on i*_nb2+j*_nb+k
    vector<Node*> _blocktrees;        // root node of search tree or null for each block, indexed on i*_nb2+j*_nb+k

    // allow our path segment generator to access our private data members
    class MySegmentGenerator;
    friend class MySegmentGenerator;

    // data members initialized during setup
    double _norm{0.};  // in the latter case, normalization factor obtained from medium component

    // data members initialized by setupSelfBefore()
    Random* _random{nullptr};
    int _numSamples{0};

    // lists of medium components of each material type;
    // list remains empty if no criteria are enabled for the corresponding material type
    vector<Medium*> _dustMedia;
    vector<Medium*> _electronMedia;
    vector<Medium*> _gasMedia;

    // flags become true if corresponding criterion is enabled
    // (i.e. configured maximum is nonzero and material type is present)
    bool _hasAny{false};
    bool _hasDustAny{false};
    bool _hasDustFraction{false};
    bool _hasDustOpticalDepth{false};
    bool _hasDustDensityDispersion{false};
    bool _hasElectronFraction{false};
    bool _hasGasFraction{false};

    // cached values for each material type (valid if corresponding flag is enabled)
    double _dustMass{0.};
    double _dustKappa{0.};
    double _electronNumber{0.};
    double _gasNumber{0.};
};

//////////////////////////////////////////////////////////////////////

#endif
