/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PLANARCUTSFORM_HPP
#define PLANARCUTSFORM_HPP

#include "GenericForm.hpp"

//////////////////////////////////////////////////////////////////////

/** PlanarCutsForm represents a generic probe form. Refer to the ProbeFormBridge class for more
    information about probes and forms.

    This particular form outputs three FITS files representing cuts through the spatial domain of
    the simulated model along three planes parallel to the coordinate planes. The user can
    configure the offset of each cut plane from the corresponding coordinate plane, the field of
    view of each cut in the relevant directions, and the number of pixels in each direction. By
    default, the cuts are centered on the origin of the model coordinate frame. If the simulation
    includes one or more media components, the default field of view is the extent of the spatial
    grid in the simulation. If there are no media in the simulation, the field of view must be
    configured explicitly.

    Each of the three FITS files contains a number of image frames corresponding to the number of
    components in a value of the quantity being probed (i.e. 1 for scalar quantities, 3 for vector
    quantities, and N for compound quantities). In case of a vector quantity, the three image
    frames representing the velocity vector components in the frame of the cut, i.e. the two
    components projected on the x and y axes of the cut and the component perpendicular to the cut,
    where positive values indicate vectors pointing away from the viewer. */
class PlanarCutsForm : public GenericForm
{
    ITEM_CONCRETE(PlanarCutsForm, GenericForm, "configurable planar cuts parallel to the coordinate planes")

        PROPERTY_DOUBLE(positionX, "the x position of the cut parallel to the yz plane")
        ATTRIBUTE_QUANTITY(positionX, "length")
        ATTRIBUTE_DEFAULT_VALUE(positionX, "0")

        PROPERTY_DOUBLE(positionY, "the y position of the cut parallel to the xz plane")
        ATTRIBUTE_QUANTITY(positionY, "length")
        ATTRIBUTE_DEFAULT_VALUE(positionY, "0")

        PROPERTY_DOUBLE(positionZ, "the z position of the cut parallel to the xy plane")
        ATTRIBUTE_QUANTITY(positionZ, "length")
        ATTRIBUTE_DEFAULT_VALUE(positionZ, "0")

        PROPERTY_DOUBLE(fieldOfViewX, "the field of view in the x direction, or zero for spatial grid bounds")
        ATTRIBUTE_QUANTITY(fieldOfViewX, "length")
        ATTRIBUTE_MIN_VALUE(fieldOfViewX, "[0")
        ATTRIBUTE_DEFAULT_VALUE(fieldOfViewX, "0")

        PROPERTY_DOUBLE(fieldOfViewY, "the field of view in the y direction, or zero for spatial grid bounds")
        ATTRIBUTE_QUANTITY(fieldOfViewY, "length")
        ATTRIBUTE_MIN_VALUE(fieldOfViewY, "[0")
        ATTRIBUTE_DEFAULT_VALUE(fieldOfViewY, "0")

        PROPERTY_DOUBLE(fieldOfViewZ, "the field of view in the z direction, or zero for spatial grid bounds")
        ATTRIBUTE_QUANTITY(fieldOfViewZ, "length")
        ATTRIBUTE_MIN_VALUE(fieldOfViewZ, "[0")
        ATTRIBUTE_DEFAULT_VALUE(fieldOfViewZ, "0")

        PROPERTY_INT(numPixelsX, "the number of pixels in the x direction")
        ATTRIBUTE_MIN_VALUE(numPixelsX, "1")
        ATTRIBUTE_MAX_VALUE(numPixelsX, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numPixelsX, "1024")

        PROPERTY_INT(numPixelsY, "the number of pixels in the y direction")
        ATTRIBUTE_MIN_VALUE(numPixelsY, "1")
        ATTRIBUTE_MAX_VALUE(numPixelsY, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numPixelsY, "1024")

        PROPERTY_INT(numPixelsZ, "the number of pixels in the z direction")
        ATTRIBUTE_MIN_VALUE(numPixelsZ, "1")
        ATTRIBUTE_MAX_VALUE(numPixelsZ, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numPixelsZ, "1024")

    ITEM_END()

public:
    /** This function causes the form to output file(s) as described in the class header for the
        quantity being probed according to the information provided by the specified
        ProbeFormBridge instance. */
    void writeQuantity(const ProbeFormBridge* bridge) const override;

    /** This static function is also used by the DefaultPlanarCutsForm. It outputs a single FITS
        file representing a cut through the spatial domain along a plane parallel to the coordinate
        plane indicated by the boolean "direction" arguments \em xd, \em yd, and \em zd, exactly
        two of which must be true. The arguments \em xp, \em yp, and \em zp specify the position of
        the cuts, the arguments \em xf, \em yf, and \em zf specify the field of view in each
        direction (or zero to indicate the spatial grid bounds), and the arguments \em xn, \em yn,
        and \em zn specify the number of pixels in each direction. */
    static void writePlanarCut(const ProbeFormBridge* bridge, bool xd, bool yd, bool zd, double xp, double yp,
                               double zp, double xf, double yf, double zf, int xn, int yn, int zn);
};

//////////////////////////////////////////////////////////////////////

#endif
