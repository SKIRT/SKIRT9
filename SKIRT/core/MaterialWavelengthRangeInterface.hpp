/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MATERIALWAVELENGTHRANGEINTERFACE_HPP
#define MATERIALWAVELENGTHRANGEINTERFACE_HPP

#include "WavelengthRangeInterface.hpp"

////////////////////////////////////////////////////////////////////

/** MaterialWavelengthRangeInterface is a pure interface. It is implemented by simulation items
    other than sources, media and instruments to indicate that the simulation item may require
    wavelength-dependent material properties for a given range of wavelengths. This allows the
    Configuration object to determine the overall wavelength range for which material properties
    must be provided in a particular simulation without requiring a compile-time dependence on the
    types of the requesting objects.

    Because the Configuration object uses this interface very early in the setup process, the
    wavelengthRange() function implementation should not rely on the receiving object being already
    set up. In other words, the wavelengthRange() function may need to call setup() on itself.

    In case the wavelengthRange() function determines at run-time (e.g. depending on
    user-configurable options) that no wavelength-dependent material properties will be requested
    by the object, the function can return a default-constructed range to indicate this.

    This interface inherits the WavelengthRangeInterface interface without extending its features,
    so it essentially serves as a wavelength interface with a special name, allowing to selectively
    search for this interface (as opposed to other interfaces inheriting WavelengthRangeInterface).
    */
class MaterialWavelengthRangeInterface : public WavelengthRangeInterface
{
protected:
    /** The empty constructor for the interface. */
    MaterialWavelengthRangeInterface() { }
};

/////////////////////////////////////////////////////////////////////////////

#endif
