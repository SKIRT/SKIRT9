/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MATERIALWAVELENGTHRANGEINTERFACE_HPP
#define MATERIALWAVELENGTHRANGEINTERFACE_HPP

#include "WavelengthRangeInterface.hpp"
class WavelengthGrid;

////////////////////////////////////////////////////////////////////

/** MaterialWavelengthRangeInterface is an interface implemented by simulation items other than
    sources, media and instruments to indicate that the simulation item may require
    wavelength-dependent material properties for a given range or set of wavelengths. This allows
    the Configuration object to determine the wavelengths for which material properties must be
    provided in a particular simulation without requiring a compile-time dependence on the types of
    the requesting objects.

    This interface extends the WavelengthRangeInterface interface, which declares the
    wavelengthRange() function, allowing to selectively search for this interface as opposed to
    other interfaces inheriting WavelengthRangeInterface.

    Because the Configuration object uses this interface very early in the setup process, the
    wavelengthRange() function implementation should not rely on the receiving object being already
    set up. In other words, the wavelengthRange() function may need to call setup() on itself or on
    the children on which its implementation depends.

    In case the wavelengthRange() function determines at run-time (e.g. depending on
    user-configurable options) that no wavelength-dependent material properties will be requested
    by the object, the function can return a default-constructed range (i.e. [0,0]) to indicate
    this.

    This interface additionally defines the materialWavelengthGrid() function, which should be
    implemented by subclasses that require material properties for the characteristic wavelengths
    of a wavelength grid. Because a default implementation is provided here for this function, this
    class is technically a mix-in class rather than a pure interface. */
class MaterialWavelengthRangeInterface : public WavelengthRangeInterface
{
protected:
    /** The empty constructor for the interface. */
    MaterialWavelengthRangeInterface() {}

public:
    /** If the receiving simulation item requires material properties at the characteristic
        wavelengths of a wavelength grid, this function returns a pointer to that wavelength grid.
        Otherwise the function returns the null pointer. The default implementation provided here
        returns the null pointer. */
    virtual WavelengthGrid* materialWavelengthGrid() const { return nullptr; }
};

/////////////////////////////////////////////////////////////////////////////

#endif
