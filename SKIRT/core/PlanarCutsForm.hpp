/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PLANARCUTSFORM_HPP
#define PLANARCUTSFORM_HPP

#include "GenericForm.hpp"
class Box;

//////////////////////////////////////////////////////////////////////

/** PlanarCutsForm represents a generic probe form. Refer to the ProbeFormBridge class for more
    information about probes and forms.

    This particular form outputs three FITS files representing cuts through the spatial domain of
    the simulated model along three planes parallel to the coordinate planes. The user can
    configure the extent of the spatial domain to be considered, the offset of each cut plane from
    the corresponding coordinate plane, and the number of pixels in each direction. By default, the
    cuts are centered on the origin of the model coordinate frame.

    Each of the three FITS files contains a number of image frames corresponding to the number of
    components in a value of the quantity being probed (i.e. 1 for scalar quantities, 3 for vector
    quantities, and N for compound quantities). In case of a vector quantity, the three image
    frames represent the velocity vector components in the frame of the cut, i.e. the two
    components projected on the x and y axes of the cut and the component perpendicular to the cut,
    where positive values indicate vectors pointing away from the viewer. */
class PlanarCutsForm : public GenericForm
{
    ITEM_CONCRETE(PlanarCutsForm, GenericForm, "configurable planar cuts parallel to the coordinate planes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(PlanarCutsForm, "Level2")

        PROPERTY_DOUBLE(minX, "the start point of the box in the X direction")
        ATTRIBUTE_QUANTITY(minX, "length")

        PROPERTY_DOUBLE(maxX, "the end point of the box in the X direction")
        ATTRIBUTE_QUANTITY(maxX, "length")

        PROPERTY_DOUBLE(minY, "the start point of the box in the Y direction")
        ATTRIBUTE_QUANTITY(minY, "length")

        PROPERTY_DOUBLE(maxY, "the end point of the box in the Y direction")
        ATTRIBUTE_QUANTITY(maxY, "length")

        PROPERTY_DOUBLE(minZ, "the start point of the box in the Z direction")
        ATTRIBUTE_QUANTITY(minZ, "length")

        PROPERTY_DOUBLE(maxZ, "the end point of the box in the Z direction")
        ATTRIBUTE_QUANTITY(maxZ, "length")

        PROPERTY_DOUBLE(positionX, "the x position of the cut parallel to the yz plane")
        ATTRIBUTE_QUANTITY(positionX, "length")
        ATTRIBUTE_DEFAULT_VALUE(positionX, "0")

        PROPERTY_DOUBLE(positionY, "the y position of the cut parallel to the xz plane")
        ATTRIBUTE_QUANTITY(positionY, "length")
        ATTRIBUTE_DEFAULT_VALUE(positionY, "0")

        PROPERTY_DOUBLE(positionZ, "the z position of the cut parallel to the xy plane")
        ATTRIBUTE_QUANTITY(positionZ, "length")
        ATTRIBUTE_DEFAULT_VALUE(positionZ, "0")

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

    /** This static function is also used by the DefaultCutsForm. It outputs a single FITS file
        representing a cut through the spatial domain along a plane parallel to the coordinate
        plane indicated by the boolean "direction" arguments \em xd, \em yd, and \em zd, exactly
        two of which must be true. The argument \em box specfies the extent of the domain to be
        considered, the arguments \em xp, \em yp, and \em zp specify the position of the cuts, and
        the arguments \em xn, \em yn, and \em zn specify the number of pixels in each direction. */
    static void writePlanarCut(const ProbeFormBridge* bridge, bool xd, bool yd, bool zd, const Box& box, double xp,
                               double yp, double zp, int xn, int yn, int zn);
};

//////////////////////////////////////////////////////////////////////

#endif
