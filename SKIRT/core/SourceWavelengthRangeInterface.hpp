/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SOURCEWAVELENGTHRANGEINTERFACE_HPP
#define SOURCEWAVELENGTHRANGEINTERFACE_HPP

#include "WavelengthRangeInterface.hpp"

////////////////////////////////////////////////////////////////////

/** SourceWavelengthRangeInterface is a pure interface. It is implemented by primary or secondary
    sources so that the source wavelength range can be obtained by other simulation items without
    requiring a compile-time dependence on the originating object's type.

    This interface inherits the WavelengthRangeInterface interface without extending its features,
    so it essentially serves as a wavelength interface with a special name, allowing to selectively
    search for this interface (as opposed to other interfaces inheriting WavelengthRangeInterface).
    */
class SourceWavelengthRangeInterface : public WavelengthRangeInterface
{
protected:
    /** The empty constructor for the interface. */
    SourceWavelengthRangeInterface() {}
};

/////////////////////////////////////////////////////////////////////////////

#endif
