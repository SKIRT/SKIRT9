/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ALLCELLSLIBRARY_HPP
#define ALLCELLSLIBRARY_HPP

#include "SpatialCellLibrary.hpp"

//////////////////////////////////////////////////////////////////////

/** The AllCellsLibrary class provides a library scheme for grouping spatial cells that has a
    separate entry for every cell, effectively disabling the library mechanism. This
    straightforward mechanism is meaningful when the corresponding calculation is fast, or as a
    reference for evaluating more complicated library mechanisms. */
class AllCellsLibrary : public SpatialCellLibrary
{
    ITEM_CONCRETE(AllCellsLibrary, SpatialCellLibrary,
                  "a library scheme that has a separate entry for every spatial cell")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked programmatically to create an AllCellsLibrary object.
        Before the constructor returns, the newly created object is hooked up as a child to the
        specified parent in the simulation hierarchy (so it will automatically be deleted), and its
        setup() function has been called. */
    explicit AllCellsLibrary(SimulationItem* parent);

    //======================== Other Functions =======================

protected:
    /** This function returns the number of entries in the library. In this class the function
        returns the number of spatial cells, i.e. \f$N_\text{entries}=N_{\text{cells}}\f$ */
    int numEntries() const override;

    /** This function returns a vector \em nv with length \f$N_{\text{cells}}\f$ that maps each
        cell index \f$m\f$ to the corresponding library entry index \f$n_m\f$. In this class the
        function returns the identity mapping. */
    vector<int> mapping(const Array& bv) const override;
};

////////////////////////////////////////////////////////////////////

#endif
