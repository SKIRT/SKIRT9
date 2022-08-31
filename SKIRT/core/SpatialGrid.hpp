/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPATIALGRID_HPP
#define SPATIALGRID_HPP

#include "Box.hpp"
#include "Position.hpp"
#include "SimulationItem.hpp"
class PathSegmentGenerator;
class Random;
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

    /** This function returns the actual or approximate diagonal of the cell with index \f$m\f$.
        For cuboidal cells, the function returns the actual diagonal. For other geometric forms, it
        returns some approximate diagonal. */
    virtual double diagonal(int m) const = 0;

    /** This function returns the index \f$m\f$ of the cell that contains the position
        \f${\bf{r}}\f$. */
    virtual int cellIndex(Position bfr) const = 0;

    /** This function returns the central location of the cell with index \f$m\f$. */
    virtual Position centralPositionInCell(int m) const = 0;

    /** This function returns a random location from the cell with index \f$m\f$. */
    virtual Position randomPositionInCell(int m) const = 0;

    /** This function creates and hands over ownership of a path segment generator (an instance of
        a PathSegmentGenerator subclass) for the spatial grid. The function must be implemented in
        each subclass to produce a path segment generator of a type corresponding to the specific
        type of spatial grid. */
    virtual std::unique_ptr<PathSegmentGenerator> createPathSegmentGenerator() const = 0;

    //================ Functions that may be implemented in subclasses ===============

public:
    /** This function outputs text data files that allow plotting the structure of the spatial
        grid. The number of data files written depends on the dimension of the spatial grid: for
        spherical symmetry only the intersection with the xy plane is written, for axial symmetry
        the intersections with the xy and xz planes are written, and for general geometries all
        three intersections are written. In the latter case, an extra file with three-dimensional
        information is written as well.

        The function is called from the SpatialGridPlotProbe class, and receives a pointer to the
        probe as its argument. The output files are called <tt>prefix_probe_grid_XXX.dat</tt>,
        where XXX is replaced by "xy", "xz", "yz" or "xyz" depending on the file under
        consideration. Within a file, each line contains two coordinates seperated by whitespace or
        is empty. Consecutive nonempty lines represent a sequence of "lineto" commands; an empty
        line marks a "moveto" command.

        The default implementation of this function invokes the write_XXX() functions as
        appropriate for the dimension of the grid. Subclasses usually override the write_XXX()
        functions, but can opt to override this function instead if that makes more sense. */
    virtual void writeGridPlotFiles(const SimulationItem* probe) const;

protected:
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
