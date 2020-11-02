/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TREESPATIALGRID_HPP
#define TREESPATIALGRID_HPP

#include "BoxSpatialGrid.hpp"
class TextOutFile;
class TreeNode;

//////////////////////////////////////////////////////////////////////

/** TreeSpatialGrid is an abstract subclass of the BoxSpatialGrid class, and represents
    three-dimensional spatial grids with cuboidal cells organized in a hierarchical tree. The
    tree's root node encloses the complete spatial domain, and nodes on subsequent levels
    recursively divide space into ever finer nodes. The depth of the tree can vary from place to
    place. The leaf nodes (those that are not further subdivided) are the actual spatial cells.

    The actual tree construction, including the choice of node type (a subclass of TreeNode) is
    delegated to each concrete subclass. This base class implements all other aspects required for
    using the grid, such as calculating paths traversing the grid. Depending on the type of
    TreeNode, the tree can become an octtree (8 children per node) or a binary tree (2 children per
    node). Other node types could be implemented, as long as they are cuboids lined up with the
    coordinate axes. */
class TreeSpatialGrid : public BoxSpatialGrid
{
    ITEM_ABSTRACT(TreeSpatialGrid, BoxSpatialGrid, "a hierarchical tree spatial grid")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** The destructor deletes all tree nodes created during setup. */
    ~TreeSpatialGrid();

protected:
    /** This function invokes the constructTree() function, to be implemented by a subclass,
        causing the tree to be constructed. The subclass returns a list of all created nodes back
        to the base class, and gives each node an identifier (ID) corresponding to its index in
        this list. Ownership of the nodes resides in the list passed back to the base class (not in
        the subclass, and not in the node hierarchy itself).

        After the subclass passes back the tree nodes, this function creates an extra vector that
        contains the node IDs of all leaf nodes, i.e. all nodes corresponding to the actual spatial
        cells. Conversely, the function also creates a vector with the cell indices of all the
        nodes, i.e. the rank \f$m\f$ of the node in the ID vector if the node is a leaf, and the
        number -1 if the node is not a leaf (and hence not a spatial cell). Finally, the function
        logs some details on the number of cells in the tree. */
    void setupSelfAfter() override;

    /** This function must be implemented in a subclass. It constructs the hierarchical tree and
        all (interconnected) nodes forming the tree. All nodes are instances of the same TreeNode
        subclass, selected by this function. The function returns a list of pointers to all created
        nodes (leaf and nonleaf). Ownership of the nodes resides in this returned list (not in the
        node hierarchy itself) and is thus handed to the caller.

        Each (leaf and nonleaf) node is given an identifier (ID) corresponding to its index in the
        list. The first node in the list is the root node of tree, i.e. the node encompasssing the
        complete spatial domain. Thus, by definition, the root node has an ID of zero.

        This function also causes the nodes to construct neighbor lists, interconnecting the nodes
        according to their spatial relationship. The neighbor lists are sorted so that the
        neighbors with the largest overlapping border area are listed first, increasing (on
        average) the probability of locating the correct neighbor early in the list. */
    virtual vector<TreeNode*> constructTree() = 0;

    //======================== Other Functions =======================

public:
    /** This function returns the number of cells in the grid. */
    int numCells() const override;

    /** This function returns the volume of the cell with index \f$m\f$. For a tree grid, it
        determines the node ID corresponding to the cell index \f$m\f$, and then simply calculates
        the volume of the corresponding cuboidal node using \f$V = \Delta x\, \Delta y\, \Delta
        z\f$. */
    double volume(int m) const override;

    /** This function returns the actuale diagonal of the cell with index \f$m\f$. For a tree grid,
        it determines the node ID corresponding to the cell index \f$m\f$, and then simply
        calculates the diagonal of the corresponding cuboidal node using \f$d = \sqrt{ (\Delta x)^2
        + (\Delta y)^2 + (\Delta z)^2 }\f$. */
    double diagonal(int m) const override;

    /** This function returns the index of the cell that contains the position \f${\bf{r}}\f$. For
        a tree grid, the search algorithm starts at the root node and selects the child node that
        contains the position. This procedure is repeated until the node is childless, i.e. until
        it is a leaf node that corresponds to an actual spatial cell. */
    int cellIndex(Position bfr) const override;

    /** This function returns the central location of the cell with index \f$m\f$. For a tree grid,
        it determines the node ID corresponding to the cell index \f$m\f$, and then calculates the
        central position in that node through \f[ \begin{split} x &= x_{\text{min}} + \frac12\,
        \Delta x \\ y &= y_{\text{min}} + \frac12\, \Delta y \\ z &= z_{\text{min}} + \frac12\,
        \Delta z \end{split} \f] */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location from the cell with index \f$m\f$. For a tree grid,
        it determines the node ID corresponding to the cell index \f$m\f$, and then calculates a
        random position in that cell through \f[ \begin{split} x &= x_{\text{min}} + {\cal{X}}_1\,
        \Delta x \\ y &= y_{\text{min}} + {\cal{X}}_2\, \Delta y \\ z &= z_{\text{min}} +
        {\cal{X}}_3\, \Delta z \end{split} \f] with \f${\cal{X}}_1\f$, \f${\cal{X}}_2\f$ and
        \f${\cal{X}}_3\f$ three uniform deviates. */
    Position randomPositionInCell(int m) const override;

    /** This function creates and hands over ownership of a path segment generator (an instance of
        a PathSegmentGenerator subclass) appropriate for a tree grid, implemented as a private
        PathSegmentGenerator subclass. The algorithm used to construct the path is described below.

        The function uses a rather straighforward algorithm. It determines the
        cell that contains the starting position, and calculates the first wall of the cell that
        will be crossed. The pathlength \f$\Delta s\f$ is determined and the current position is
        moved to a new position along this path, a tiny fraction further than \f$\Delta s\f$, \f[
        \begin{split} x_{\text{new}} &= x_{\text{current}} + (\Delta s + \epsilon)\,k_x \\
        y_{\text{new}} &= y_{\text{current}} + (\Delta s + \epsilon)\,k_y \\ z_{\text{new}} &=
        z_{\text{current}} + (\Delta s + \epsilon)\,k_z \end{split} \f] where \f[ \epsilon =
        10^{-12} \sqrt{x_{\text{max}}^2 + y_{\text{max}}^2 + z_{\text{max}}^2} \f] By adding this
        small extra bit, we ensure that the new position is now within the next cell, and we can
        repeat this exercise. This loop is terminated when the next position is outside the grid.

        To determine the cell index of the "next cell" in this algorithm, the function uses the
        neighbor lists constructed for each tree node during setup. */
    std::unique_ptr<PathSegmentGenerator> createPathSegmentGenerator() const override;

    /** This function writes the topology of the tree to the specified text file in a simple,
        proprietary format. After a brief descriptive header, it writes lines that each contain
        just a single integer number. The first line specifies the number of children for each
        nonleaf node (2 for a binary tree, 8 for an octtree, or 0 if the root node has not been
        subdivided). The second line contains 1 if the root node is subdivided, or 0 if not. The
        following lines similarly contain 1 or 0 indicating subdivision for any children of the
        preceding node, recursively, in a depth-first traversal of the tree. */
    void writeTopology(TextOutFile* outfile) const;

protected:
    /** This function writes the intersection of the grid with the xy plane to the specified
        SpatialGridPlotFile object. */
    void write_xy(SpatialGridPlotFile* outfile) const override;

    /** This function writes the intersection of the grid with the xz plane to the specified
        SpatialGridPlotFile object. */
    void write_xz(SpatialGridPlotFile* outfile) const override;

    /** This function writes the intersection of the grid with the yz plane to the specified
        SpatialGridPlotFile object. */
    void write_yz(SpatialGridPlotFile* outfile) const override;

    /** This function writes 3D information for the cells up to a certain level in the grid
        structure to the specified SpatialGridPlotFile object. The output is limited to some
        predefined number of cells to keep the line density in the output plot within reason. */
    void write_xyz(SpatialGridPlotFile* outfile) const override;

private:
    /** This function returns a pointer to the root node of the tree. */
    TreeNode* root() const;

    /** This function returns a pointer to the node corresponding to cell index \f$m\f$. It just
        reads the node ID of the \f$m\f$'th leaf cell from the precalculated ID vector and returns
        the corresponding pointer of the tree vector. */
    TreeNode* nodeForCellIndex(int m) const;

    /** This function returns the cell index \f$m\f$ of a node in the tree. It just obtains the
        node ID from the node and determines the corresponding cell index from the precalculated
        cell index vector. */
    int cellIndexForNode(const TreeNode* node) const;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    double _eps{0.};           // a small fraction relative to the spatial extent of the grid
    vector<TreeNode*> _nodev;  // list of all nodes in the tree; holds ownership; first item is root node
                               // node id in each node corresponds to index in this vector
    vector<int> _cellindexv;   // cell index m corresponding to each node in nodev; -1 for nonleaf nodes
    vector<int> _idv;          // node id (or equivalently, index in nodev) for each cell (i.e. leaf node)

    // allow our path segment generator to access our private data members
    class MySegmentGenerator;
    friend class MySegmentGenerator;
};

//////////////////////////////////////////////////////////////////////

#endif
