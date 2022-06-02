/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef VORONOIMESHSNAPSHOT_HPP
#define VORONOIMESHSNAPSHOT_HPP

#include "Array.hpp"
#include "Snapshot.hpp"
class PathSegmentGenerator;
class SiteListInterface;
class SpatialGridPath;

////////////////////////////////////////////////////////////////////

/** A VoronoiMeshSnapshot object represents a Voronoi tessellation or \em mesh of a cuboidal
    spatial domain (given a list of generating sites) and offers several facilities related to this
    mesh. A Voronoi mesh partitions the domain in convex polyhedra. Consider a given set of points
    in the domain, called \em sites. For each site there will be a corresponding region consisting
    of all points in the domain closer to that site than to any other. These regions are called
    Voronoi \em cells, and together they form the Voronoi tessellation of the domain.

    The VoronoiMeshSnapshot class serves two key purposes: (1) represent snapshot data produced by
    a hydrodynamical simulation and imported from a column text file, defining a primary source or
    a transfer medium distribution, and (2) implement an unstructured spatial grid that discretizes
    the spatial domain as a basis for the radiative transfer simulation itself.

    To support the first use case (imported snapshot data), the VoronoiMeshSnapshot class is based
    on the Snapshot class; it uses the facilities offered there to configure and help read the
    snapshot data, and it implements all functions in the general Snapshot public interface. In
    addition it offers functionality that is specific to this snapshot type, such as, for example,
    the requirement to configure the spatial extent of the domain. In this use case, the client
    employs the default constructor and then proceeds to configure the snapshot object as described
    in the Snapshot class header.

    To support the second use case (spatial grid for radiative transfer), the class offers
    specialty constructors that accept a list of the generating sites from various sources,
    including a programmatically generated list or a user-supplied text column file. Note that this
    input file includes just the site coordinates, while a snapshot data file would include
    additional properties for each site. The specialty constructors automatically complete the
    configuration sequence of the object, so that the getters can be used immediately after
    construction.

    Once an VoronoiMeshSnapshot object has been constructed and fully configured, its data members
    are no longer modified. Consequently all getters are re-entrant.

    Using the Voro++ library
    ------------------------

    To build the Voronoi tessellation, the buildMesh() function in this class uses the code in the
    \c voro subfolder of the SKIRT code hierarchy, which is taken from the Voro++ library written
    by Chris H. Rycroft (Harvard University/Lawrence Berkeley Laboratory) at
    https://github.com/chr1shr/voro (git commit 122531f) with minimal changes to avoid compiler
    warnings.

    A distinguishing feature of the Voro++ library is that it carries out cell-based calculations,
    computing the Voronoi cell for each site individually, rather than computing the Voronoi
    tessellation as a global network of vertices and edges. It is therefore particularly
    well-suited for applications that require cell-based properties such as the cell volume, the
    centroid position or the number of faces or vertices. Equally important in the context of
    SKIRT, it is easy to distribute the work over parallel execution threads because, after setting
    up a common search data structure holding all sites, the calculations for the various cells are
    mutually independent.

    The Voro++ approach also has an important drawback. Because cells are handled independently of
    each other, floating pointing rounding errors sometimes cause inconsistencies where the results
    for one cell do not properly match the corresponding results for a neighboring cell. This most
    often happens when the generating sites are very close to each other and/or form certain
    hard-to-calculate geometries. These problems manifest themselves either as empty cells or as
    asymmetries in the neighbor lists. To handle these situations, the buildMesh() function removes
    some sites from the input list as described below.

    Firstly, before actually building the Voronoi tessellation, and in addition to removing sites
    that lie outside of the spatial domain, the function discards sites that lie closer to any
    previously listed site than \f$10^{-12}\f$ times the diagonal of the spatial domain. This
    preventative measure essentially removes any "degenerate" sites from the input list. Given the
    extremely small distances, it is very unlikely that the physcis of the input model would be
    affected by this operation.

    Secondly, after all Voronoi cells have been calculated, the function verifies the cell
    properties. If any of the cells have a zero volume or if any of the cell neighbors are not
    mutual, the involved sites are discarded and the Voronoi construction starts anew from scratch
    with the reduced list of sites. After a maximum of 5 attempts, the function throws a fatal
    error. In practice, this will hopefully never happen. Discarding incorrectly calculated cells
    perhaps incurs a slightly higher risk of changing the physcis of the input model. However, in
    practice it seems that these issues mostly occur in regions of high site density, so that the
    errors should be fairly limited. */
class VoronoiMeshSnapshot : public Snapshot
{
    //================= Construction - Destruction =================

public:
    /** The default constructor initializes the snapshot in an invalid state; see the description
        of the required calling sequence in the Snapshot class header. */
    VoronoiMeshSnapshot();

    /** The destructor releases any data structures allocated by this class. */
    ~VoronoiMeshSnapshot();

    //========== Reading ==========

public:
    /** This function reads the snapshot data from the input file, honoring the options set through
        the configuration functions, stores the data for later use, and closes the file by calling
        the base class Snapshot::readAndClose() function. Sites located outside of the domain and
        sites that are too close to another site are discarded. Sites with an associated
        temperature above the cutoff temperature (if one has been configured) are assigned a
        density value of zero, so that the corresponding cell has zero mass (regardless of the
        imported mass/density properties).

        The function calls the private buildMesh() function to build the Voronoi mesh based on the
        imported site positions. If the snapshot configuration requires the ability to determine
        the density at a given spatial position, the function also calls the private buildSearch()
        function to create a data structure that accelerates locating the cell containing a given
        point.

        During its operation, the function logs some statistical information about the imported
        snapshot and the resulting data structures. */
    void readAndClose() override;

    //========== Configuration ==========

public:
    /** This function sets the extent of the spatial domain for the Voronoi mesh. When using the
        default constructor, this function must be called during configuration. There is no
        default; failing to set the extent of the domain results in undefined behavior. */
    void setExtent(const Box& extent);

    /** This function configures the snapshot to skip construction of the actual Voronoi
        tessellation and instead use a search tree across all sites. It should be called only if
        (1) the snapshot has been configured to import both a mass/number density column \em and a
        volume-integrated mass/number column, and (2) the snapshot will not be required to generate
        random positions or trace paths. Violating these conditions will result in undefined
        behavior. */
    void foregoVoronoiMesh();

    //========== Specialty constructors ==========

public:
    /** This constructor reads the site positions from the specified text column file. The input
        file must contain three columns specifying the x,y,z coordinates. The default unit is
        parsec, which can be overridden by providing column header info in the file. The
        constructor completes the configuration for the object (but without importing mass density
        information or setting a mass density policy) and calls the private buildMesh() and
        buildSearch() functions to create the relevant data structures.

        The \em item argument specifies a simulation item in the hierarchy of the caller (usually
        the caller itself) used to retrieve context such as an appropriate logger. The \em extent
        argument specifies the extent of the domain as a box lined up with the coordinate axes.
        Sites located outside of the domain and sites that are too close to another site are
        discarded. The \em filename argument specifies the name of the input file, including
        filename extension but excluding path and simulation prefix. If the \em relax argument is
        true, the function performs a single relaxation step on the site positions. */
    VoronoiMeshSnapshot(const SimulationItem* item, const Box& extent, string filename, bool relax);

    /** This constructor obtains the site positions from a SiteListInterface instance. The
        constructor completes the configuration for the object (but without importing mass density
        information or setting a mass density policy) and calls the private buildMesh() and
        buildSearch() functions to create the relevant data structures.

        The \em item argument specifies a simulation item in the hierarchy of the caller (usually
        the caller itself) used to retrieve context such as an appropriate logger. The \em extent
        argument specifies the extent of the domain as a box lined up with the coordinate axes.
        Sites located outside of the domain and sites that are too close to another site are
        discarded. The \em sli argument specifies an object that provides the SiteListInterface
        interface from which to obtain the site positions. If the \em relax argument is true, the
        function performs a single relaxation step on the site positions. */
    VoronoiMeshSnapshot(const SimulationItem* item, const Box& extent, SiteListInterface* sli, bool relax);

    /** This constructor obtains the site positions from a programmatically prepared list. The
        constructor completes the configuration for the object (but without importing mass density
        information or setting a mass density policy) and calls the private buildMesh() and
        buildSearch() functions to create the relevant data structures.

        The \em item argument specifies a simulation item in the hierarchy of the caller (usually
        the caller itself) used to retrieve context such as an appropriate logger. The \em extent
        argument specifies the extent of the domain as a box lined up with the coordinate axes.
        Sites located outside of the domain and sites that are too close to another site are
        discarded. The \em sites argument specifies the list of site positions. If the \em relax
        argument is true, the function performs a single relaxation step on the site positions. */
    VoronoiMeshSnapshot(const SimulationItem* item, const Box& extent, const vector<Vec>& sites, bool relax);

    //=========== Private construction ==========

private:
    /** Private class to hold the information about a Voronoi cell that is relevant for calculating
        paths and densities; see the buildMesh() function. */
    class Cell;

    /** Private class to hold a node in the internal binary search tree; see the buildTree()
        and buildSearch() functions. */
    class Node;

    /** Given a list of generating sites (represented as partially initialized Cell
        objects), this private function builds the Voronoi tessellation and stores the
        corresponding cell information, including any properties relevant for supporting the
        interrogation capabilities offered by this class. All other data (such as Voronoi cell
        vertices, edges and faces) are discarded. In practice, the function adds the sites to a
        Voro++ container, computes the Voronoi cells one by one, and copies the relevant cell
        information (such as the list of neighboring cells) from the Voro++ data structures into
        its own.

        If the \em relax argument is true, the function performs a single relaxation step on the
        site positions using Lloyd's algorithm (Lloyd 1982; Du, Faber and Gunzburger 1999, SIAM
        review 41.4, pp 637-676; Dobbels 2017, master thesis). An intermediate Voronoi tessellation
        is built using the original site positions, and subsequently each site position is replaced
        by the centroid (mass center) of the corresponding cell. The final tessellation is then
        constructed with these adjusted site positions, which are distributed more uniformly,
        thereby avoiding overly elongated cells in the Voronoi tessellation. Relaxation can be
        quite time-consuming because the Voronoi tessellation must be constructed twice. */
    void buildMesh(bool relax);

    /** This private function calculates the volumes for all cells without using the Voronoi mesh.
        It assumes that both mass and mass density columns are being imported. */
    void calculateVolume();

    /** This private function calculates the densities and (cumulative) masses for all cells, and
        logs some statistics. The function assumes that the cell volumes have been calculated,
        either by building a Voronoi tessellation, or by deriving the volume from mass and mass
        density columns being imported. */
    void calculateDensityAndMass();

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

    /** This private function builds a data structure that allows accelerating the operation of the
        cellIndex() function without using the Voronoi mesh. The domain is not partitioned in
        blocks. The function builds a single binary search tree on all cell sites (see for example
        <a href="http://en.wikipedia.org/wiki/Kd-tree">en.wikipedia.org/wiki/Kd-tree</a>). */
    void buildSearchSingle();

    /** This private function returns true if the given point is closer to the site with index m
        than to the sites with indices ids. */
    bool isPointClosestTo(Vec r, int m, const vector<int>& ids) const;

    //====================== Output =====================

public:
    /** This function outputs grid plot files as described for the SpatialGridPlotProbe. The
        function reconstructs the Voronoi tesselation in order to produce the coordinates of the
        Voronoi cell vertices. */
    void writeGridPlotFiles(const SimulationItem* probe) const;

    //=========== Interrogation ==========

public:
    /** This function returns the extent of the spatial domain as configured through the
        setExtent() function. */
    Box extent() const override;

    /** This function returns the number of sites (or, equivalently, cells) in the snapshot. */
    int numEntities() const override;

    /** This function returns the position of the site with index \em m. If the index is out of
        range, the behavior is undefined. */
    Position position(int m) const override;

    /** This function returns the centroid of the Voronoi cell with index \em m. If the index is
        out of range, the behavior is undefined. */
    Position centroidPosition(int m) const;

    /** This function returns the volume of the Voronoi cell with index \em m. If the index is out
        of range, the behavior is undefined. */
    double volume(int m) const override;

    /** This function returns the bounding box (enclosing cuboid lined up with the coordinate axes)
        of the Voronoi cell with index \em m. If the index is out of range, the behavior is
        undefined. */
    Box extent(int m) const;

    /** This function returns the mass density associated with the cell with index \em m. If no
        density policy has been set or no mass information is being imported, or if the index is
        out of range, the behavior is undefined. */
    double density(int m) const override;

    /** This function returns the mass density represented by the snapshot at a given point
        \f${\bf{r}}\f$, or equivalently, the mass density associated with the cell containing the
        given point. If the point is outside the domain, the function returns zero. If no density
        policy has been set or no mass information is being imported, or if the search data
        structures used by the cellIndex() function were not created during construction, the
        behavior is undefined. */
    double density(Position bfr) const override;

    /** This function returns the total mass represented by the snapshot, in other words the sum of
        the masses of all cells. If no density policy has been set or no mass information is being
        imported, the behavior is undefined. */
    double mass() const override;

    /** This function returns a random position drawn uniformly from the (polyhedron) volume of the
        cell with index \em m. If the index is out of range, the behavior is undefined.

        The function generates uniformly distributed random points in the enclosing cuboid until
        one happens to be inside the cell. The candidate point is inside the cell if it is closer
        to the cell's site position than to any neighbor cell's site positions. */
    Position generatePosition(int m) const override;

    /** This function returns a random position within the spatial domain of the snapshot, drawn
        from the mass density distribution represented by the snapshot. The function first selects
        a random cell from the discrete probability distribution formed by the respective cell
        masses, and then generates a random position uniformly from the volume of that cell. If no
        density policy has been set or no mass information is being imported, the behavior is
        undefined. */
    Position generatePosition() const override;

    /** This function returns the cell index \f$0\le m \le N_{cells}-1\f$ for the cell containing
        the specified point \f${\bf{r}}\f$. If the point is outside the domain, the function
        returns -1. By definition of a Voronoi tesselation, the closest site position determines
        the Voronoi cell containing the specified point.

        The function uses the search data structures created by the private BuildSearch() function
        to accelerate its operation. It computes the appropriate block index from the coordinates
        of the specified point, which provides a list of Voronoi cells possibly overlapping the
        point. If there is a search tree for this block, the function uses it to locate the nearest
        point in \f${\cal{O}}(\log N)\f$ time. Otherwise it calculates the distance from the
        specified point to the site positions for each of the possibly overlapping cells,
        determining the nearest one in linear time. For a small number of cells this is more
        efficient than using the search tree.

        If the search data structures were not created during construction (which happens when
        using the default constructor without configuring a mass density policy), invoking the
        cellIndex() function causes undefined behavior. */
    int cellIndex(Position bfr) const;

protected:
    /** This function returns a reference to an array containing the imported properties (in column
        order) for the cell with index \f$0\le m \le N_\mathrm{ent}-1\f$. If the index is out of
        range, the behavior is undefined. */
    const Array& properties(int m) const override;

    /** This function returns the index \f$0\le m \le N_\mathrm{ent}-1\f$ of the cell containing
        the specified point \f${\bf{r}}\f$, or -1 if the point is outside the domain, if there
        are no cells in the snapshot, or if the search data structures were not created. */
    int nearestEntity(Position bfr) const override;

public:
    /** This function sets the specified entity collection to the cell containing the specified
        point \f${\bf{r}}\f$, or to the empty collection if the point is outside the domain or if
        there are no cells in the snapshot. If the search data structures were not created,
        invoking this function causes undefined behavior. */
    void getEntities(EntityCollection& entities, Position bfr) const override;

    /** This function replaces the contents of the specified entity collection by the set of cells
        crossed by the specified path with starting point \f${\bf{r}}\f$ and direction
        \f${\bf{k}}\f$. The weight of a cell is given by the length of the path segment inside the
        cell. If the path does not cross the spatial domain of the snapshot, the collection will be
        empty. If the search data structures were not created, invoking this function causes
        undefined behavior. */
    void getEntities(EntityCollection& entities, Position bfr, Direction bfk) const override;

    //====================== Path construction =====================

public:
    /** This function creates and hands over ownership of a path segment generator appropriate for
        the adaptive mesh spatial grid, implemented as a private PathSegmentGenerator subclass. The
        algorithm used to construct the path is described below.

        In the first stage, the function checks whether the start point is inside the domain. If
        so, the current point is simply initialized to the start point. If not, the function
        computes the path segment to the first intersection with one of the domain walls and moves
        the current point inside the domain. Finally the function determines the current cell, i.e.
        the cell containing the current point.

        In the second stage, the function loops over the algorithm that computes the exit point
        from the current cell, i.e. the intersection of the ray formed by the current point and
        the path direction with the current cell's boundary. By the nature of Voronoi cells, this
        algorithm also produces the ID of the neigboring cell without extra cost. If an exit point
        is found, the loop adds a segment to the output path, updates the current point and the
        current cell, and continues to the next iteration. If the exit is towards a domain wall,
        the path is complete and the loop is terminated. If no exit point is found, which shouldn't
        happen too often, this must be due to computational inaccuracies. In that case, no path
        segment is added, the current point is advanced by a small amount, and the new current cell
        is determined by calling the function whichcell().

        The algorithm that computes the exit point has the following input data:
        <TABLE>
        <TR><TD>\f$\bf{r}\f$</TD>               <TD>%Position of the current point</TD></TR>
        <TR><TD>\f$\bf{k}\f$</TD>               <TD>%Direction of the ray, normalized</TD></TR>
        <TR><TD>\f$m_r\f$</TD>                  <TD>ID of the cell containing the current point</TD></TR>
        <TR><TD>\f$m_i,\;i=0\ldots n-1\f$</TD>  <TD>IDs of the neighboring cells (\f$m_i>=0\f$)
                                                    or domain walls (\f$m_i<0\f$)</TD></TR>
        <TR><TD>\f$\mathbf{p}(m)\f$</TD>        <TD>Site position for cell \f$m\f$ (implicit)</TD></TR>
        <TR><TD>\f$x_\text{min},x_\text{max},y_\text{min},y_\text{max},z_\text{min},z_\text{max}\f$</TD>
                                                <TD>Domain boundaries (implicit)</TD></TR>
        </TABLE>
        where the domain wall IDs are defined as follows:
        <TABLE>
        <TR><TD><B>Domain wall ID</B></TD>      <TD><B>Domain wall equation</B></TD></TR>
        <TR><TD>-1</TD>                         <TD>\f$x=x_\text{min}\f$</TD></TR>
        <TR><TD>-2</TD>                         <TD>\f$x=x_\text{max}\f$</TD></TR>
        <TR><TD>-3</TD>                         <TD>\f$y=y_\text{min}\f$</TD></TR>
        <TR><TD>-4</TD>                         <TD>\f$y=y_\text{max}\f$</TD></TR>
        <TR><TD>-5</TD>                         <TD>\f$z=z_\text{min}\f$</TD></TR>
        <TR><TD>-6</TD>                         <TD>\f$z=z_\text{max}\f$</TD></TR>
        </TABLE>

        The line containing the ray can be written as \f$\mathcal{L}(s)=\mathbf{r}+s\,\mathbf{k}\f$
        with \f$s\in\mathbb{R}\f$. The exit point can similarly be written as
        \f$\mathbf{q}=\mathbf{r}+s_q\,\mathbf{k}\f$ with \f$s_q>0\f$, and the distance covered within
        the cell is given by \f$s_q\f$. The ID of the cell next to the exit point is denoted \f$m_q\f$
        and is easily determined as explained below.

        The algorithm that computes the exit point proceeds as follows:
         - Calculate the set of values \f$\{s_i\}\f$ for the intersection points between the line
           \f$\mathcal{L}(s)\f$ and the planes defined by the neighbors or walls \f$m_i\f$
           (see below for details on the intersection calculation).
         - Select the smallest nonnegative value \f$s_q=\min\{s_i|s_i>0\}\f$ in the set to determine
           the exit point. The neighbor or wall ID with the same index \f$i\f$ determines the
           corresponding \f$m_q\f$.
         - If there is no nonnegative value in the set, no exit point has been found.

        To calculate \f$s_i\f$ for a regular neighbor \f$m_i\geq 0\f$, intersect the
        line \f$\mathcal{L}(s)=\mathbf{r}+s\,\mathbf{k}\f$ with the plane bisecting the points
        \f$\mathbf{p}(m_i)\f$ and \f$\mathbf{p}(m_r)\f$. An unnormalized vector perpendicular to
        this plane is given by \f[\mathbf{n}=\mathbf{p}(m_i)-\mathbf{p}(m_r)\f]
        and a point on the plane is given by
        \f[\mathbf{p}=\frac{\mathbf{p}(m_i)+\mathbf{p}(m_r)}{2}.\f] The equation
        of the plane can then be written as \f[\mathbf{n}\cdot(\mathbf{x}-\mathbf{p})=0.\f]
        Substituting \f$\mathbf{x}=\mathbf{r}+s_i\,\mathbf{k}\f$ and solving for \f$s_i\f$ provides
        \f[s_i=\frac{\mathbf{n}\cdot(\mathbf{p}-\mathbf{r})}{\mathbf{n}\cdot\mathbf{k}}.\f]
        If \f$\mathbf{n}\cdot\mathbf{k}=0\f$ the line and the plane are parallel and there
        is no intersection. In that case no \f$s_i\f$ is added to the set of candidate
        exit points.

        To calculate \f$s_i\f$ for a wall \f$m_i<0\f$, substitute the appropriate normal and
        position vectors for the wall plane in this last formula. For example, for the left wall
        with \f$m_i=-1\f$ one has \f$\mathbf{n}=(-1,0,0)\f$ and \f$\mathbf{p}=(x_\text{min},0,0)\f$
        so that \f[s_i=\frac{x_\text{min}-r_x}{k_x}.\f]
    */
    std::unique_ptr<PathSegmentGenerator> createPathSegmentGenerator() const;

    //======================== Data Members ========================

private:
    // data members initialized during configuration
    Box _extent;                     // the spatial domain of the mesh
    double _eps{0.};                 // small fraction of extent
    bool _foregoVoronoiMesh{false};  // true if using search tree instead of Voronoi tessellation

    // data members initialized when processing snapshot input and further completed by BuildMesh()
    vector<Cell*> _cells;  // cell objects, indexed on m

    // data members initialized when processing snapshot input, but only if a density policy has been set
    Array _rhov;       // density for each cell (not normalized)
    Array _cumrhov;    // normalized cumulative density distribution for cells
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
};

////////////////////////////////////////////////////////////////////

#endif
