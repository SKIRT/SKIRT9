/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef COLORDECORATOR_HPP
#define COLORDECORATOR_HPP

#include "ShapeDecorator.hpp"

////////////////////////////////////////////////////////////////////

/** The ColorDecorator class represents a decorator that sets the color for the shape or group of
    shapes being decorated. */
class ColorDecorator : public ShapeDecorator
{
    /** This enumeration lists the colors supported by the ColorDecorator. Select 'Custom' to allow
        specifying arbitrary red, green, and blue intensities. */
    ENUM_DEF(Color, White, Red, Green, Blue, Cyan, Magenta, Yellow, Black, Custom)
        ENUM_VAL(Color, White, "the white color")
        ENUM_VAL(Color, Red, "the red color")
        ENUM_VAL(Color, Green, "the green color")
        ENUM_VAL(Color, Blue, "the blue color")
        ENUM_VAL(Color, Cyan, "the cyan color")
        ENUM_VAL(Color, Magenta, "the magenta color")
        ENUM_VAL(Color, Yellow, "the yellow color")
        ENUM_VAL(Color, Black, "the black color")
        ENUM_VAL(Color, Custom, "a custom RGB color")
    ENUM_END()

    ITEM_CONCRETE(ColorDecorator, ShapeDecorator, "a decorator that sets the color of a shape or group of shapes")
        PROPERTY_ENUM(color, Color, "the color of the decorated shape or group of shapes")
        PROPERTY_DOUBLE_LIST(rgb, "the red, green, blue intensities of the decorated shape or group of shapes")
        ATTRIBUTE_MIN_VALUE(rgb, "0")
        ATTRIBUTE_MAX_VALUE(rgb, "1")
        ATTRIBUTE_DEFAULT_VALUE(rgb, "0.5, 0.5, 0.5")
        ATTRIBUTE_RELEVANT_IF(rgb, "colorCustom")
    ITEM_END()

    // ================== Drawing ==================

protected:
    /** This function sets the line width in the current graphics state by calling the setColor()
        function. */
    void paintSelf() override;
};

////////////////////////////////////////////////////////////////////

#endif
