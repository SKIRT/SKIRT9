/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ADAPTIVEMESHSNAPSHOT_HPP
#define ADAPTIVEMESHSNAPSHOT_HPP

#include "Array.hpp"
#include "Snapshot.hpp"
class PathSegmentGenerator;
class SpatialGridPath;

////////////////////////////////////////////////////////////////////

/** An AdaptiveMeshSnapshot object represents a three-dimensional Adaptive Mesh Refinement (AMR)
    grid and offers several related interrogation facilities. As implemented by this class, an
    adaptive mesh recursively partitions a cuboidal spatial domain into cuboidal (sub)cells. In
    contrast to an octree, which has a fixed 2x2x2 subdivision scheme, an adaptive mesh uses an
    arbitrary scheme that can (and must) be specified for each subdivision seperately. As a result,
    an adaptive mesh can provide high resolution (i.e. small cells) in areas where it matters,
    without consuming memory in areas where the resolution can be lower, in a very flexible manner.

    The primary objective of the AdaptiveMeshSnapshot class is to represent snapshot data produced
    by a hydrodynamical simulation and imported from a column text file, defining a primary source
    or a transfer medium distribution. To support this use case, the class is based on the Snapshot
    class; it uses the facilities offered there to configure and help read the snapshot data, and
    it implements all functions in the general Snapshot public interface. In addition it offers
    functionality that is specific to this snapshot type, such as, for example, the requirement to
    configure the spatial extent of the domain. A client should employ the default
    AdaptiveMeshSnapshot constructor and configure the snapshot object as described in the Snapshot
    class header.

    Additionally, once a snapshot has been imported, the adaptive mesh defined by the snapshot can
    in some cases also be used to discretize the spatial domain as a basis for the radiative
    transfer simulation itself. To support this use case, the class offers the capability to trace
    a linear path through the adaptive mesh.

    Once an AdaptiveMeshSnapshot object has been constructed and fully configured, its data is no
    longer modified. Consequently all getters are re-entrant.

    Tree structure, Morton ordering, and file format
    ------------------------------------------------

    The adaptive mesh data structure is organized in a tree. Each tree node represents a cubodial
    portion of the domain, called its extent. The root node's extent is the complete domain. A
    nonleaf node distributes its extent over its child nodes using a regular linear grid. The
    number of subdivisions is defined separately for each node and can differ for each spatial
    dimension. A leaf node represents a cell in which physical quantities are considered to be
    constant. Collectively the leaf nodes form a partition of the domain, i.e. their extents cover
    the complete domain without overlapping one another.

    The cells in this three-dimensional data structure can be arranged in a linear sequence using
    Morton ordering. This ordering is obtained by performing a depth-first traversal of the tree,
    where each nonleaf node outputs its children in the order x-first, then y, then z.

    When reading an adaptive mesh snapshot from a text column input file, each line in the file
    describes a particular tree node (nonleaf or leaf), and the lines are given in Morton order.
    Specifically, each line in the file can be one of the following types:

     - Comment: lines with a crosshatch (#) as the first non-whitespace character, lines
       containing only whitespace and empty lines are ignored (and do not count in the Morton order).
     - Nonleaf: a nonleaf line has an exclamation mark (!) as the first non-whitespace character,
       followed by optional whitespace and then three whitespace-separated positive integer numbers
       \f$N_x,N_y,N_z\f$. These three numbers specify the number of child nodes carried by this node
       in each spatial direction. The child nodes are on a regular linear grid as described above.
     - Leaf: a leaf node contains one or more whitespace-separated floating point numbers reflecting
       the physical quantities associated with the leaf node, depending on the user configuration
       settings for the snapshot. The default units for these quantities can be overridden by
       including columns header info in the file as described in the TextInFile class header.

    Note that the input file does not include the physical extent of the domain; this information
    must be provided as part of the user configuration for the snapshot.

    Below is an example of an adaptive mesh import file. For illustrative purposes:

      - the mesh is assumed to have a single cell in the z direction;
      - the value in the first colum indicates the Morton order index of the leaf cell;
      - the values in the remaining two columns indicate the x and y coordinates of
        the cell's center using a domain of size 4. x 3. with one corner at the origin.

    \verbatim
    # Example adaptive mesh import data file
    #
    ! 4 3 1
     0 0.50 0.50
     1 1.50 0.50
     2 2.50 0.50
     3 3.50 0.50
     4 0.50 1.50
     5 1.50 1.50
    ! 2 2 1
     6 2.25 1.25
     7 2.75 1.25
     8 2.25 1.75
    ! 2 2 1
     9 2.63 1.63
    10 2.88 1.63
    11 2.63 1.88
    12 2.88 1.88
    ! 2 2 1
    13 3.25 1.25
    14 3.75 1.25
    15 3.25 1.75
    16 3.75 1.75
    17 0.50 2.50
    18 1.50 2.50
    ! 2 2 1
    19 2.25 2.25
    20 2.75 2.25
    21 2.25 2.75
    22 2.75 2.75
    ! 2 2 1
    23 3.25 3.25
    24 3.75 3.25
    25 3.25 3.75
    26 3.75 3.75
    \endverbatim

*/
class AdaptiveMeshSnapshot : public Snapshot
{
    //================= Construction - Destruction =================

public:
    /** The default constructor initializes the snapshot in an invalid state; see the description
        of the required calling sequence in the Snapshot class header. */
    AdaptiveMeshSnapshot();

    /** The destructor releases any data structures allocated by this class. */
    ~AdaptiveMeshSnapshot();

    //========== Reading ==========

public:
    /** This function reads the snapshot data from the input file, honoring the options set through
        the configuration functions, stores the data for later use, and closes the file by calling
        the base class Snapshot::readAndClose() function.

        Cells with an associated temperature above the cutoff temperature (if one has been
        configured) are assigned a density value of zero, so that the cell has zero mass regardless
        of the imported mass/density properties.

        The function logs some statistical information about the imported snapshot and the
        resulting data structures. */
    void readAndClose() override;

    //========== Configuration ==========

public:
    /** This function sets the extent of the spatial domain for the adaptive mesh. It must be
        called during configuration, before the readAndClose() function is invoked. There is no
        default extent; failing to set the extent of the domain results in undefined behavior. */
    void setExtent(const Box& extent);

    /** This function adds neighbor information to all leaf nodes in the adaptive mesh. If should
        be called after the readAndClose() function has completed its operation, and before the
        createPathSegmentGenerator() function is invoked. Specifically, the function causes each
        leaf node to remember its most likely neighbor at each of its six walls. This information,
        while optional, substantially accelerates path construction. */
    void addNeighbors();

    //=========== Interrogation ==========

public:
    /** This function returns the extent of the spatial domain as configured through the
        setExtent() function. */
    Box extent() const override;

    /** This function returns the number of leaf cells in the adaptive mesh snapshot. */
    int numEntities() const override;

    /** This function returns the position of the center of the leaf cell with index \em m. If the
        index is out of range, the behavior is undefined. */
    Position position(int m) const override;

    /** This function returns the volume of the leaf cell with index \em m. If the index is out of
        range, the behavior is undefined. */
    double volume(int m) const override;

    /** This function returns the diagonal of the leaf cell with index \em m. If the index is out of
        range, the behavior is undefined. */
    double diagonal(int m) const;

    /** This function returns the bounding box (which is by definition lined up with the coordinate
        axes) of the leaf cell with index \em m. If the index is out of range, the behavior is
        undefined. */
    Box extent(int m) const;

    /** This function returns the mass density associated with the leaf cell with index \em m. If
        no density policy has been set or no mass information is being imported, or if the index is
        out of range, the behavior is undefined. */
    double density(int m) const override;

    /** This function returns the mass density represented by the snapshot at a given point
        \f${\bf{r}}\f$, or equivalently, the mass density associated with the leaf cell containing
        the given point. If the point is outside the domain, the function returns zero. If no
        density policy has been set or no mass information is being imported, the behavior is
        undefined. */
    double density(Position bfr) const override;

    /** This function returns the total mass represented by the snapshot, in other words the sum of
        the masses of all leaf cells. If no density policy has been set or no mass information is
        being imported, the behavior is undefined. */
    double mass() const override;

    /** This function returns a random position drawn uniformly from the (cuboidal) extent of the
        leaf cell with index \em m. If the index is out of range, the behavior is undefined. */
    Position generatePosition(int m) const override;

    /** This function returns a random position within the spatial domain of the snapshot, drawn
        from the mass density distribution represented by the snapshot. The function first selects
        a random leaf cell from the discrete probability distribution formed by the respective cell
        masses, and then generates a random position uniformly from the volume of that cell. If no
        density policy has been set or no mass information is being imported, the behavior is
        undefined. */
    Position generatePosition() const override;

    /** This function returns the leaf cell index \f$0\le m \le N_{cells}-1\f$ for the cell
        containing the specified point \f${\bf{r}}\f$. If the point is outside the domain, the
        function returns -1. The function recursively searches the adaptive mesh tree until it
        finds the appropriate leaf cell. */
    int cellIndex(Position bfr) const;

protected:
    /** This function returns a reference to an array containing the imported properties (in column
        order) for the cell with index \f$0\le m \le N_\mathrm{ent}-1\f$. If the index is out of
        range, the behavior is undefined. */
    const Array& properties(int m) const override;

    /** This function returns the index \f$0\le m \le N_\mathrm{ent}-1\f$ of the cell containing
        the specified point \f${\bf{r}}\f$, or -1 if the point is outside the domain. */
    int nearestEntity(Position bfr) const override;

public:
    /** This function sets the specified entity collection to the cell containing the specified
        point \f${\bf{r}}\f$, or to the empty collection if the point is outside the domain. */
    void getEntities(EntityCollection& entities, Position bfr) const override;

    /** This function replaces the contents of the specified entity collection by the set of cells
        crossed by the specified path with starting point \f${\bf{r}}\f$ and direction
        \f${\bf{k}}\f$. The weight of a cell is given by the length of the path segment inside the
        cell. If the path does not cross the spatial domain of the snapshot, the collection will be
        empty. */
    void getEntities(EntityCollection& entities, Position bfr, Direction bfk) const override;

    //====================== Path construction =====================

public:
    /** This function creates and hands over ownership of a path segment generator appropriate for
        the adaptive mesh spatial grid, implemented as a private PathSegmentGenerator subclass.

        The algorithm used to construct the path is fairly straightforward because all cells are
        cuboids lined up with the coordinate axes. The information added by the addNeighbors()
        function significantly accelerates path construction. */
    std::unique_ptr<PathSegmentGenerator> createPathSegmentGenerator() const;

    //======================== Data Members ========================

private:
    // data members initialized during configuration
    Box _extent;      // the spatial domain of the mesh
    double _eps{0.};  // small fraction of extent

    // data members initialized when processing snapshot input
    class Node;
    Node* _root{nullptr};  // root node representing the complete domain
    vector<Node*> _cells;  // leaf nodes indexed on m

    // data members initialized when processing snapshot input, but only if a density policy has been set
    Array _rhov;       // density for each cell (not normalized)
    Array _cumrhov;    // normalized cumulative density distribution for cells
    double _mass{0.};  // total effective mass

    // allow our path segment generator to access our private data members
    class MySegmentGenerator;
    friend class MySegmentGenerator;
};

////////////////////////////////////////////////////////////////////

#endif
