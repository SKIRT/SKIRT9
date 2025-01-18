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

/** A CylindricalCellSnapshot object represents the data in a snapshot imported from a column text
    file and discretized using a cylindrical coordinate-grid.

    Each line in the text file represents a cylindrical cell lined up with the cylindrical
    coordinate axes. The columns specify the coordinates of the bordering planes and cylinders (see
    the CylCell class documentation for more information) along with properties such as the density
    in the cell.

    The intention is for the cells to define a partition of the spatial domain, and usually they
    will, however this is not enforced. When two or more cells overlap at a given position, the
    properties for that position will be taken from the cell that is listed first in the imported
    file. When no cells overlap a given position, the density at that position is considered to be
    zero. To avoid thin slices of zero density between cells, the coordinates for common walls or
    corners in neighboring cells should be identical.

    This class is based on the Snapshot class; it uses the facilities offered there to configure
    and read the snapshot data, and it implements all functions in the general Snapshot public
    interface. If the snapshot configuration requires the ability to determine the density at a
    given spatial position, an effort is made to accelerate the search for the cell containing that
    position even for a large number of cells. */
class CylindricalCellSnapshot : public Snapshot
{
    //================= Construction - Destruction =================

    //========== Reading ==========

public:
    /** This function reads the snapshot data from the input file, honoring the options set through
        the configuration functions, stores the data for later use, and closes the file by calling
        the base class Snapshot::close() function.

        Cells with an associated temperature above the cutoff temperature (if one has been
        configured) are assigned a density value of zero, so that they have zero mass regardless of
        the imported mass/density properties. Note that we cannot simply ignore these cells because
        an empty cell may overlap and thus hide (a portion of) a nonempty cell later in the list.

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
