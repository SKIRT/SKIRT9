/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef WAVELENGTHRANGEINTERFACE_HPP
#define WAVELENGTHRANGEINTERFACE_HPP

#include "Range.hpp"

////////////////////////////////////////////////////////////////////

/** WavelengthRangeInterface is an abstract base class for pure interfaces that can be implemented
    by simulation items wishing to communicate a wavelength range to other simulation items without
    requiring a compile-time dependence on the originating object's type. */
class WavelengthRangeInterface
{
protected:
    /** The empty constructor for the interface. */
    WavelengthRangeInterface() {}

public:
    /** The empty destructor for the interface. */
    virtual ~WavelengthRangeInterface() {}

    /** This function returns the wavelength range represented by the object implementing the interface. */
    virtual Range wavelengthRange() const = 0;
};

/////////////////////////////////////////////////////////////////////////////

#endif
