/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef WAVELENGTHRANGEINTERFACE_HPP
#define WAVELENGTHRANGEINTERFACE_HPP

#include "Range.hpp"

////////////////////////////////////////////////////////////////////

/** WavelengthRangeInterface is a pure interface. It is implemented by the primary source system so
    that the primary source wavelength range for the simulation can be obtained by other simulation
    items without requiring a dependence on the source system class. */
class WavelengthRangeInterface
{
protected:
    /** The empty constructor for the interface. */
    WavelengthRangeInterface() { }

public:
    /** The empty destructor for the interface. */
    virtual ~WavelengthRangeInterface() { }

    /** This function returns the wavelength range represented by the object implementing the interface. */
    virtual Range wavelengthRange() const = 0;
};

/////////////////////////////////////////////////////////////////////////////

#endif
