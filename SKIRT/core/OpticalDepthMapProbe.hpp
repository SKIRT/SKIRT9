/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef OPTICALDEPTHMAPPROBE_HPP
#define OPTICALDEPTHMAPPROBE_HPP

#include "AbstractWavelengthProbe.hpp"
#include "AllSkyProjection.hpp"

////////////////////////////////////////////////////////////////////

/** OpticalDepthMapProbe outputs a FITS file for each material type with an all-sky optical depth
    map at the given wavelength as seen from the given position towards the given direction. Each
    map has a 2:1 aspect ratio (with square pixels) and uses the user-selected projection to
    project the complete sky on the rectangular image. Pixels outside of the mapped region are set
    to zero. The files are named <tt>prefix_probe_MM_tau.fits</tt> where MM is replaced by the
    material type. */
class OpticalDepthMapProbe : public AbstractWavelengthProbe
{
    /** The enumeration type indicating when probing occurs. */
    ENUM_DEF(ProbeAfter, Setup, Run)
        ENUM_VAL(ProbeAfter, Setup, "after setup")
        ENUM_VAL(ProbeAfter, Run, "after the complete simulation run")
    ENUM_END()

    ITEM_CONCRETE(OpticalDepthMapProbe, AbstractWavelengthProbe,
                  "all-sky optical depth map as seen from given location")
        ATTRIBUTE_TYPE_DISPLAYED_IF(OpticalDepthMapProbe, "Level2&Medium&SpatialGrid")

        PROPERTY_ITEM(projection, AllSkyProjection, "the projection used for mapping the sky to a rectangle")
        ATTRIBUTE_DEFAULT_VALUE(projection, "HammerAitoffProjection")

        PROPERTY_INT(numPixelsY, "the number of image pixels in the vertical (shortest) direction")
        ATTRIBUTE_MIN_VALUE(numPixelsY, "25")
        ATTRIBUTE_MAX_VALUE(numPixelsY, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numPixelsY, "250")

        PROPERTY_DOUBLE(observerX, "the position of the observer, x component")
        ATTRIBUTE_QUANTITY(observerX, "length")
        ATTRIBUTE_DEFAULT_VALUE(observerX, "0")

        PROPERTY_DOUBLE(observerY, "the position of the observer, y component")
        ATTRIBUTE_QUANTITY(observerY, "length")
        ATTRIBUTE_DEFAULT_VALUE(observerY, "0")

        PROPERTY_DOUBLE(observerZ, "the position of the observer, z component")
        ATTRIBUTE_QUANTITY(observerZ, "length")
        ATTRIBUTE_DEFAULT_VALUE(observerZ, "0")

        PROPERTY_DOUBLE(crossX, "the position of the crosshair, x component")
        ATTRIBUTE_QUANTITY(crossX, "length")
        ATTRIBUTE_DEFAULT_VALUE(crossX, "1")

        PROPERTY_DOUBLE(crossY, "the position of the crosshair, y component")
        ATTRIBUTE_QUANTITY(crossY, "length")
        ATTRIBUTE_DEFAULT_VALUE(crossY, "0")

        PROPERTY_DOUBLE(crossZ, "the position of the crosshair, z component")
        ATTRIBUTE_QUANTITY(crossZ, "length")
        ATTRIBUTE_DEFAULT_VALUE(crossZ, "0")

        PROPERTY_DOUBLE(upX, "the upwards direction, x component")
        ATTRIBUTE_QUANTITY(upX, "length")
        ATTRIBUTE_DEFAULT_VALUE(upX, "0")

        PROPERTY_DOUBLE(upY, "the upwards direction, y component")
        ATTRIBUTE_QUANTITY(upY, "length")
        ATTRIBUTE_DEFAULT_VALUE(upY, "0")

        PROPERTY_DOUBLE(upZ, "the upwards direction, z component")
        ATTRIBUTE_QUANTITY(upZ, "length")
        ATTRIBUTE_DEFAULT_VALUE(upZ, "1")

        PROPERTY_ENUM(probeAfter, ProbeAfter, "when to probe the medium state")
        ATTRIBUTE_DEFAULT_VALUE(probeAfter, "Setup")
        ATTRIBUTE_DISPLAYED_IF(probeAfter, "DynamicState")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that all attribute values have been appropriately set. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function performs probing after setup. It produces output only if the \em
        probeAfter property is set to Setup. */
    void probeSetup() override;

    /** This function performs probing after all photon packets have been emitted and detected. It
        produces output only if the \em probeAfter property is set to Run. */
    void probeRun() override;

private:
    /** This function performs the probing; it is called from probeSetup() or probeRun() depending
        on the value of the \em probeAfter property. */
    void probe();

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation
    const double& _Ox{_observerX};
    const double& _Oy{_observerY};
    const double& _Oz{_observerZ};
    const double& _Cx{_crossX};
    const double& _Cy{_crossY};
    const double& _Cz{_crossZ};
    const double& _Ux{_upX};
    const double& _Uy{_upY};
    const double& _Uz{_upZ};
};

////////////////////////////////////////////////////////////////////

#endif
