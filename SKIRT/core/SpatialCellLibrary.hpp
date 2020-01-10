/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPATIALCELLLIBRARY_HPP
#define SPATIALCELLLIBRARY_HPP

#include "Array.hpp"
#include "SimulationItem.hpp"

//////////////////////////////////////////////////////////////////////

/** SpatialCellLibrary is an abstract base class for providing spatial cell grouping schemes called
    libraries. A library maps a number of spatial cells to a single library entry based on some
    predefined scheme, such as similarities in the stored radiation field or other cell properties.
    This in turn allows a client to perform common (possibly approximate) calculations for all of
    the cells mapped to the same library entry, trading accuracy for speed. This base class offers
    just an interface that must be implemented in each subclass.

    The SpatialCellLibrary class and its subclasses support the implementation of a library
    mechanism as described in Baes et al. (2011, ApJS, 196, 22). Instead of calculating the
    emission spectrum individually for every spatial cell in the system, a library is constructed,
    a template spectrum is calculated for each library entry, and these templates are used for all
    spatial cells mapped to the corresponding library entry. Obviously, the spectral templates in
    the library should be chosen/constructed in such a way that they can represent the whole range
    of actual spectra encountered in the simulation. In other words, the library should span the
    entire parameter space of interstellar radiation fields. Different subclasses of the
    SpatialCellLibrary class achieve this goal to different degrees of sophistication, with a
    better coverage of the parameter space typically at the cost of a more resource-intensive
    library construction. */
class SpatialCellLibrary : public SimulationItem
{
    ITEM_ABSTRACT(SpatialCellLibrary, SimulationItem, "a library scheme for grouping spatial cells")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the number of entries \f$N_\text{entries}\f$ in the library. Library
        entries might remain unused in the sense that none of the indices returned by the mapping()
        function point to these entries. This function must be implemented by each subclass. */
    virtual int numEntries() const = 0;

    /** This function returns a vector \em nv with length \f$N_\text{cells}\f$ that maps each cell
        index \f$m\f$ to the corresponding library entry index \f$n_m\f$. The indices in the
        returned vector are in the range \f$[-1,N_\text{entries}-1]\f$. An index value of -1
        indicates that the cell is not included in the mapping and should not be used (for example,
        because the cell will produce a negligible amount of emission or no emission at all).

        The argument array \em bv with length \f$N_\text{cells}\f$ provides a value passed from the
        caller for each cell in the spatial grid. If the value is zero, the cell will not be used
        by the caller regardless of the returned mapping index, and thus it can safely be omitted
        from the mapping (i.e. given an index of -1). If the \em bv value for the cell is nonzero,
        the caller plans to use the cell, but it will still refrain from doing so if the library
        decides not to map the cell (i.e. give it an index of -1).

        This function must be implemented by each subclass. */
    virtual vector<int> mapping(const Array& bv) const = 0;
};

////////////////////////////////////////////////////////////////////

#endif
