/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPATIALGRID_HPP
#define SPATIALGRID_HPP

#include "SimulationItem.hpp"
#include "Box.hpp"
#include "Position.hpp"
class Random;
class SpatialGridPath;
class SpatialGridPlotFile;

//////////////////////////////////////////////////////////////////////

/** The SpatialGrid class is an abstract base class for grids that tessellate the spatial domain of
    the simulation. Each position in the computational domain corresponds to a single spatial cell.
    A SpatialGrid subclass instance represents only purely geometric properties, i.e.\ it contains
    no information on the actual distribution of material over the grid. */
class SpatialGrid : public SimulationItem
{
    ITEM_ABSTRACT(SpatialGrid, SimulationItem, "a spatial grid")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function caches the simulation's random generator for use by subclasses. */
    void setupSelfBefore() override;

    //================ Functions that must be implemented in subclasses ===============

public:
    /** This function returns the dimension of the grid, which depends on its (lack of) symmetry. A
        value of 1 means spherical symmetry, 2 means axial symmetry and 3 means none of these
        symmetries. */
    virtual int dimension() const = 0;

    /** This function returns the number of cells in the grid. */
    virtual int numCells() const = 0;

    /** This function returns the bounding box that encloses the grid. */
    virtual Box boundingBox() const = 0;

    /** This function returns the volume of the cell with index \f$m\f$. */
    virtual double volume(int m) const = 0;

    /** This function returns the index \f$m\f$ of the cell that contains the position
        \f${\bf{r}}\f$. */
    virtual int cellIndex(Position bfr) const = 0;

    /** This function returns the central location of the cell with index \f$m\f$. */
    virtual Position centralPositionInCell(int m) const = 0;

    /** This function returns a random location from the cell with index \f$m\f$. */
    virtual Position randomPositionInCell(int m) const = 0;

    /** This function calculates a path through the grid. The SpatialGridPath object passed as an
        argument specifies the starting position \f${\bf{r}}\f$ and the direction \f${\bf{k}}\f$
        for the path. The data on the calculated path are added back into the same object. This
        consists of three vectors: the first one lists the cell indices \f$m\f$ of all the cells
        crossed by the path, the second lists the path length \f$\Delta s\f$ covered in each of
        these cells, and the third lists the accumulated path length \f$s\f$ until the end of each
        cell is encountered. */
    virtual void path(SpatialGridPath* path) const = 0;

    //================ Functions that may be implemented in subclasses ===============

public:
    /** This function writes the intersection of the grid with the xy plane to the specified
        SpatialGridPlotFile object. The default implementation does nothing. */
    virtual void write_xy(SpatialGridPlotFile* outfile) const;

    /** This function writes the intersection of the grid with the xz plane to the specified
        SpatialGridPlotFile object. The default implementation does nothing. */
    virtual void write_xz(SpatialGridPlotFile* outfile) const;

    /** This function writes the intersection of the grid with the yz plane to the specified
        SpatialGridPlotFile object. The default implementation does nothing. */
    virtual void write_yz(SpatialGridPlotFile* outfile) const;

    /** This function writes 3D information for all or part of the cells in the grid structure to
        the specified SpatialGridPlotFile object. The default implementation does nothing. */
    virtual void write_xyz(SpatialGridPlotFile* outfile) const;

    //======================== Other Functions =======================

protected:
    /** This function returns the simulation's random generator as a service to subclasses. */
    Random* random() const { return _random; }

    //======================== Data Members ========================

private:
    // data member initialized during setup
    Random* _random{nullptr};
};

//////////////////////////////////////////////////////////////////////

#endif
