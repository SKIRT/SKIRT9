/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CARTESIANCELLINTERFACE_HPP
#define CARTESIANCELLINTERFACE_HPP

#include "Box.hpp"

////////////////////////////////////////////////////////////////////

/** CartesianCellInterface is a pure interface implemented by spatial grids that have Cartesian
    cells, or in other words, cells that are identical to their axis-aligned bounding box. The
    interface offers a way to obtain the bounding box representing a given cell. */
class CartesianCellInterface
{
protected:
    /** The empty constructor for the interface. */
    CartesianCellInterface() {}

public:
    /** The empty destructor for the interface. */
    virtual ~CartesianCellInterface() {}

    /** This function returns the number of cells in the grid. */
    virtual int numCells() const = 0;

    /** This function returns the box defining the cell with given zero-based index. If the index
        is out of range, the behavior is undefined. */
    virtual Box cellBox(int m) const = 0;
};

/////////////////////////////////////////////////////////////////////////////

#endif
