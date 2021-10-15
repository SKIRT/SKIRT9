/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CELLSNAPSHOT_HPP
#define CELLSNAPSHOT_HPP

#include "Snapshot.hpp"

////////////////////////////////////////////////////////////////////

/** A CellSnapshot object represents the data in a Cartesian mesh-based snapshot imported from a
    column text file. Each line in the text file represents a cuboidal cell lined up with the
    coordinate axes, specifying the coordinates of the lower-left and upper-right corners of the
    cell along with properties such as the density in the cell. The intention is for the cells to
    define a partition of the spatial domain, and usually they will, however this is not enforced.
    When two or more cells overlap at a given position, the properties for that position will be
    taken from the cell that is listed first in the imported file. When no cells overlap a given
    position, the density at that position is considered to be zero. To avoid thin slices of zero
    density between cells, the coordinates for common walls or corners in neighboring cells should
    be identical.

    This class is an alternative to the AdaptiveMeshSnapshot class, which requires listing cells in
    Morton order. This can be hard to accomplish, especially when extracting a particular subdomain
    from a larger mesh. In contrast, the CellSnapshot class allows cells to be listed in arbitrary
    order and allows the spatial partition to be imperfect. For example, the outer border of the
    domain can be ragged, with some cells extending farther out than others.

    This class is based on the Snapshot class; it uses the facilities offered there to configure
    and read the snapshot data, and it implements all functions in the general Snapshot public
    interface. If the snapshot configuration requires the ability to determine the density at a
    given spatial position, an effort is made to accelerate the search for the cell containing that
    position even for a large number of cells. */
class CellSnapshot : public Snapshot
{
    //================= Construction - Destruction =================

public:
    /** The destructor deletes the search data structure, if it was constructed. */
    ~CellSnapshot();

    //========== Reading ==========

public:
    /** This function reads the snapshot data from the input file, honoring the options set through
        the configuration functions, stores the data for later use, and closes the file by calling
        the base class Snapshot::readAndClose() function.

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

    /** This function returns the position of the center of the cell with index \em m. If the index
        is out of range, the behavior is undefined. */
    Position position(int m) const override;

    /** This function returns the metallicity of the cell with index \em m. If the metallicity is
        not being imported, or the index is out of range, the behavior is undefined. */
    double metallicity(int m) const override;

    /** This function returns the metallicity of the cell containing the specified point
        \f${\bf{r}}\f$. If the point is not inside any cell, the function returns -1. If the
        metallicity is not being imported, the behavior is undefined. */
    double metallicity(Position bfr) const override;

    /** This function returns the temperature of the cell with index \em m. If the temperature is
        not being imported, or the index is out of range, the behavior is undefined. */
    double temperature(int m) const override;

    /** This function returns the temperature of the cell containing the specified point
        \f${\bf{r}}\f$. If the point is not inside any cell, the function returns -1. If the
        temperature is not being imported, the behavior is undefined. */
    double temperature(Position bfr) const override;

    /** This function returns the velocity of the cell with index \em m. If the velocity is not
        being imported, or the index is out of range, the behavior is undefined. */
    Vec velocity(int m) const override;

    /** This function returns the velocity of the cell containing the specified point
        \f${\bf{r}}\f$. If the point is not inside any cell, the function returns zero velocity. If
        the velocity is not being imported, the behavior is undefined. */
    Vec velocity(Position bfr) const override;

    /** This function returns the velocity dispersion of the cell with index \f$0\le m \le
        N_\mathrm{ent}-1\f$. If the velocity dispersion is not being imported, or the index is out
        of range, the behavior is undefined. */
    double velocityDispersion(int m) const override;

    /** This function returns the velocity dispersion of the cell containing the specified point
        \f${\bf{r}}\f$. If the point is not inside any cell, the function returns zero dispersion.
        If the velocity dispersion is not being imported, the behavior is undefined. */
    double velocityDispersion(Position bfr) const override;

    /** This function returns the magnetic field vector of the cell with index \em m. If the
        magnetic field is not being imported, or the index is out of range, the behavior is
        undefined. */
    Vec magneticField(int m) const override;

    /** This function returns the magnetic field vector of the cell containing the specified point
        \f${\bf{r}}\f$. If the point is not inside any cell, the function returns a zero magnetic
        field. If the magnetic field is not being imported, the behavior is undefined. */
    Vec magneticField(Position bfr) const override;

    /** This function stores the parameters of the cell with index \em m into the given array. If
        parameters are not being imported, or the index is out of range, the behavior is undefined.
        */
    void parameters(int m, Array& params) const override;

    /** This function stores the parameters of the cell containing the specified point
        \f${\bf{r}}\f$ into the given array. If the point is not inside any cell, the function
        returns the appropriate number of zero parameter values. If parameters are not being
        imported, the behavior is undefined. */
    void parameters(Position bfr, Array& params) const override;

    /** This function returns the mass density of the cell containing the specified point
        \f${\bf{r}}\f$. If the point is not inside any cell, the function returns zero. If no
        density policy has been set or no mass information is being imported, the behavior is
        undefined. */
    double density(Position bfr) const override;

    /** This function returns the total mass represented by the snapshot, in other words the sum of
        the masses of all cells. If no density policy has been set or no mass information is being
        imported, the behavior is undefined. */
    double mass() const override;

    /** This function returns a random position drawn uniformly from the cell with index \em m. If
        the index is out of range, the behavior is undefined. */
    Position generatePosition(int m) const override;

    /** This function returns a random position within the spatial domain of the snapshot, drawn
        from the mass density distribution represented by the snapshot. The function first selects
        a random cell from the discrete probability distribution formed by the respective cell
        masses, and then generates a random position within that cell. If no density policy has
        been set or no mass information is being imported, the behavior is undefined. */
    Position generatePosition() const override;

    //======================== Data Members ========================

private:
    // data members initialized when reading the input file
    vector<Array> _propv;  // cell properties as imported

    // data members initialized after reading the input file if a density policy has been set
    Array _rhov;       // density for each cell (not normalized)
    Array _cumrhov;    // normalized cumulative density distribution for cells
    double _mass{0.};  // total effective mass

    // data members initialized after reading the input file if a density policy has been set
    class CellGrid;
    CellGrid* _grid{nullptr};  // smart grid for locating the cell at a given location
};

////////////////////////////////////////////////////////////////////

#endif
