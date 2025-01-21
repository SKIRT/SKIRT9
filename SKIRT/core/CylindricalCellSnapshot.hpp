/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CYLINDRICALCELLSNAPSHOT_HPP
#define CYLINDRICALCELLSNAPSHOT_HPP

#include "BoxSearch.hpp"
#include "CylCell.hpp"
#include "Snapshot.hpp"

////////////////////////////////////////////////////////////////////

/** A CylindricalCellSnapshot object represents the data in a 2D or 3D snapshot imported from a
    column text file and discretized using cylindrical coordinates.

    <em>3D data</em>

    By default, each line in the text file represents a 3D cylindrical cell lined up with the
    cylindrical coordinate axes. The columns specify the coordinates of the bordering planes and
    cylinders, along with properties such as density. Specifically, each cell is bordered by:

    - two vertical cylinders centered on the origin defined by radii \f$0 \le R_\text{min} \le
    R_\text{max}\f$,

    - two meridional half-planes (with \f$R>0\f$) defined by azimuth angles \f$-\pi \le
    \varphi_\text{min} \le \varphi_\text{max} \le \pi\f$ with \f$\varphi_\text{max} -
    \varphi_\text{min} \le \pi\f$,

    - two horizontal planes defined by \f$z_\text{min} \le z_\text{max}\f$.

    \note Because of the limitations on the range of \f$\varphi\f$, a Cylcell cannot straddle the
    negative x-axis of the Cartesian model coordinate system, and it cannot span more than half of
    the azimuth circle.

    The intention is for the cells to define a partition of the spatial domain, and usually they
    will, however this is not enforced. When two or more cells overlap at a given position, the
    properties for that position will be taken from the cell that is listed first in the imported
    file. When no cells overlap a given position, the density at that position is considered to be
    zero. To avoid thin slices of zero density between cells, the coordinates for common walls or
    corners in neighboring cells should be identical.

    <em>Velocities</em>

    If velocity or magnetic field vectors are being imported, the class converts these from
    cylindrical to Cartesian coordinates using \f[\begin{aligned} v_\text{x} &= v_{R}\cos\varphi -
    v_\varphi\sin\varphi \\ v_\text{y} &= v_{R}\sin\varphi + v_\varphi\cos\varphi \\ v_\text{z} &=
    v_\text{z} \end{aligned}\f] where \f$\varphi = (\varphi_\text{min}+\varphi_\text{max})/2\f$ is
    the central angle of the corresponding cell.

    <em>2D data</em>

    If the auto-revolve feature is enabled, each line in the text file represents a 2D cell in the
    meridional plane with \f$\varphi=0\f$. The cells are defined using just the radial and vertical
    borders \f$0 \le R_\text{min} \le R_\text{max}\f$ and \f$z_\text{min} \le z_\text{max}\f$.
    After reading the text file, these 2D cells will automatically be revolved around the z-axis
    using a user-specified number of \f$\varphi\f$ bins.

    To enable the auto-revolve feature, the number of auto-revolve bins must be set to at least 2.
    Also, all \f$\varphi_\text{min}\f$ and \f$\varphi_\text{max}\f$ values in the input file must
    be exactly zero. A nonzero value in these columns will trigger a fatal error.

    \note It is \em not allowed to omit the \f$\varphi_\text{min}, \varphi_\text{max}\f$ columns.

    If the 2D input file specifies an integrated mass type (as opposed to a mass density), the mass
    of a 2D cell is evenly distributed over the revolved 3D bins.

    For a model that includes velocities, the velocities will be converted from cylindrical to
    Cartesian coordinates as described above after revolving to a 3D representation. Because the
    result of this conversion depends on the azimuth angle, the user must set a sufficiently high
    number of auto-revolve bins to properly resolve the revolved velocity field.

    <em>Implementation</em>

    This class is based on the Snapshot class; it uses the facilities offered there to configure
    and read the snapshot data, and it implements all functions in the general Snapshot public
    interface. If the snapshot configuration requires the ability to determine the density at a
    given spatial position, an effort is made to accelerate the search for the cell containing that
    position even for a large number of cells. */
class CylindricalCellSnapshot : public Snapshot
{
    //================= Construction - Destruction =================

public:
    /** This function sets the number of auto-revolve azimuth bins, enabling the auto-revolve
        feature if the number of bins is at least 2 and all azimuth angles in the imported data are
        zero (see the class header for more information). */
    void setNumAutoRevolveBins(int numBins);

    //========== Reading ==========

public:
    /** This function reads the snapshot data from the input file, honoring the options set through
        the configuration functions, stores the data for later use, and closes the file by calling
        the base class Snapshot::close() function.

        Cells with an associated temperature above the cutoff temperature (if one has been
        configured) are assigned a density value of zero, so that they have zero mass regardless of
        the imported mass/density properties. Note that we cannot simply ignore these cells because
        an empty cell may overlap and thus hide (a portion of) a nonempty cell later in the list.

        If velocity and/or magnetic field vectors are being imported, the function converts these
        from cylindrical to Cartesian coordinates as described in the class header.

        The function also logs some statistical information about the import. If the snapshot
        configuration requires the ability to determine the density at a given spatial position,
        this function builds a data structure that accelerates the search for the appropriate cell.
        */
    void readAndClose() override;

    //=========== Interrogation ==========

public:
    /** This function returns the bounding box lined up with the coordinate axes surrounding all
        cells. */
    Box extent() const override;

    /** This function returns the number of cells in the snapshot. */
    int numEntities() const override;

    /** This function returns the volume of the cell with index \em m. If the index is out of
        range, the behavior is undefined. */
    double volume(int m) const override;

    /** This function returns the mass density associated with the cell with index \em m. If no
        density policy has been set or no mass information is being imported, or if the index is
        out of range, the behavior is undefined. */
    double density(int m) const override;

    /** This function returns the mass density of the cell containing the specified point
        \f${\bf{r}}\f$. If the point is not inside any cell, the function returns zero. If no
        density policy has been set or no mass information is being imported, the behavior is
        undefined. */
    double density(Position bfr) const override;

    /** This function returns the total mass represented by the snapshot, in other words the sum of
        the masses of all cells. If no density policy has been set or no mass information is being
        imported, the behavior is undefined. */
    double mass() const override;

    /** This function returns the position of the center of the cell with index \em m. If the index
        is out of range, the behavior is undefined. */
    Position position(int m) const override;

    /** This function returns a random position drawn uniformly from the cell with index \em m. If
        the index is out of range, the behavior is undefined. */
    Position generatePosition(int m) const override;

    /** This function returns a random position within the spatial domain of the snapshot, drawn
        from the mass density distribution represented by the snapshot. The function first selects
        a random cell from the discrete probability distribution formed by the respective cell
        masses, and then generates a random position within that cell. If no density policy has
        been set or no mass information is being imported, the behavior is undefined. */
    Position generatePosition() const override;

protected:
    /** This function returns a reference to an array containing the imported properties (in column
        order) for the cell with index \f$0\le m \le N_\mathrm{ent}-1\f$. If the index is out of
        range, the behavior is undefined. */
    const Array& properties(int m) const override;

public:
    /** This function sets the specified entity collection to the cell containing the specified
        point \f${\bf{r}}\f$, or to the empty collection if the point is outside the domain or if
        there are no cells in the snapshot. If the search data structures were not created,
        invoking this function causes undefined behavior. */
    void getEntities(EntityCollection& entities, Position bfr) const override;

    /** This function replaces the contents of the specified entity collection by the set of cells
        that overlap the specified path with starting point \f${\bf{r}}\f$ and direction
        \f${\bf{k}}\f$. The weight of a cell is given by the length of the path segment inside the
        cell. If the path does not overlap any cells, the collection will be empty. If the search
        data structures were not created, invoking this function causes undefined behavior. */
    void getEntities(EntityCollection& entities, Position bfr, Direction bfk) const override;

    //======================== Data Members ========================

private:
    // data members initialized during configuration
    int _numAutoRevolveBins{0};  // must be 2 or more to enable auto-revolve feature

    // data members initialized when reading the input file
    vector<Array> _propv;    // cell properties as imported
    vector<CylCell> _cellv;  // imported coordinates converted to a cylindrical cell object

    // data members initialized after reading the input file if a density policy has been set
    Array _rhov;        // density for each cell (not normalized)
    Array _cumrhov;     // normalized cumulative density distribution for cells
    double _mass{0.};   // total effective mass
    BoxSearch _search;  // search structure for locating cells
};

////////////////////////////////////////////////////////////////////

#endif
