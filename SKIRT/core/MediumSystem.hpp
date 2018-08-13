/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEDIUMSYSTEM_HPP
#define MEDIUMSYSTEM_HPP

#include "SimulationItem.hpp"
#include "Array.hpp"
#include "Medium.hpp"
#include "Table.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the MediumSystem class represents a complete medium system, which is the
    superposition of one or more transfer media. Each individual medium represents a spatial
    density distribution and defines the material properties of the medium at each location. While
    the specific material properties may vary with location, the fundamental material type must be
    the same throughout the spatial domain for each medium.

    In addition to the media input model, the MediumSystem class includes the spatial grid that
    tessellates the spatial domain of the simulation into cells, and manages the medium state for
    each spatial cell in this grid.

    TODO: add more info on managing the medium state. */
class MediumSystem : public SimulationItem
{
    ITEM_CONCRETE(MediumSystem, SimulationItem, "a medium system")

    PROPERTY_ITEM_LIST(media, Medium, "the transfer media")
        ATTRIBUTE_DEFAULT_VALUE(media, "GeometricMedium")
        ATTRIBUTE_OPTIONAL(media)

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function XXX. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the medium system, which depends on the (lack of)
        symmetry in the geometries of its components. A value of 1 means spherical symmetry, 2
        means axial symmetry and 3 means none of these symmetries. The medium with the least
        symmetry (i.e. the highest dimension) determines the result for the whole system. */
    int dimension() const;

    /** This function returns the number of media in the medium system. */
    int numMedia() const;

    /** This function returns the number of cells in the spatial grid. */
    int numCells() const;

    //======================== Data Members ========================

private:
    // initialized during setup
    Array _volumev;     // volume of each cell (indexed on m)
    Table<2> _nvv;      // number density for each cell and each dust component (indexed on m,h)
};

////////////////////////////////////////////////////////////////

#endif
