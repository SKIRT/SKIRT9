/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef WAVELENGTHRANGEINTERFACE_HPP
#define WAVELENGTHRANGEINTERFACE_HPP

#include "Range.hpp"

////////////////////////////////////////////////////////////////////

/** WavelengthRangeInterface is a pure interface. It is implemented by primary or secondary sources
    so that the source wavelength range can be obtained by other simulation items without requiring
    a dependence on the originating object. */
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
