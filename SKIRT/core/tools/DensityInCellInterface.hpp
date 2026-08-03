/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DENSITYINCELLINTERFACE_HPP
#define DENSITYINCELLINTERFACE_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** DensityInCellInterface is a pure interface. It is implemented by spatial grids that offer a
    fast way to obtain the density of a grid cell, given the cell index \em m, for each medium
    component \em h. If the MediumSystem class detects a spatial grid that implements this
    interface, it will call the numberDensity() function in this interface rather than sampling the
    density distribution in a number of random points. In some special cases, for example when the
    grid is lined up with some cell structure in the medium geometry, this can dramatically enhance
    performance. */
class DensityInCellInterface
{
protected:
    /** The empty constructor for the interface. */
    DensityInCellInterface() {}

public:
    /** The empty destructor for the interface. */
    virtual ~DensityInCellInterface() {}

    /** This function returns the number density for medium component \em h in the spatial grid
        cell with index \em m. */
    virtual double numberDensity(int h, int m) const = 0;
};

/////////////////////////////////////////////////////////////////////////////

#endif
