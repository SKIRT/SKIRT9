/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SHAPECANVAS_HPP
#define SHAPECANVAS_HPP

#include "Shape.hpp"

////////////////////////////////////////////////////////////////////

/** The ShapeCanvas class is a shape item that manages a canvas on which the shapes in a hierarchy
    can be drawn, and provides facilities to perform such drawing. The root of a run-time shape
    hierarchy must be an instance of the ShapeCanvas class. */
class ShapeCanvas : public ShapeItem
{
    ITEM_CONCRETE(ShapeCanvas, ShapeItem, "the root of a shape hierarchy, managing the canvas")
        PROPERTY_STRING(savePath, "the path of the output file when saving this canvas")
        ATTRIBUTE_REQUIRED_IF(savePath, "false")
        ATTRIBUTE_DISPLAYED_IF(savePath, "false")
        PROPERTY_ITEM(shape, Shape, "the top-level shape in this hierarchy")
        ATTRIBUTE_DEFAULT_VALUE(shape, "ShapeGroup")
    ITEM_END()

    // ================== Drawing ==================

public:
    /** This function paints the shape hierarchy held by the receiving item onto its canvas, and
        then saves the result to the specified output file. Call the paintAndSave() function of a
        shape canvas rather than the paint() function inherited from the Shape class. */
    void paintAndSave(string outPath);

protected:
    /** This function pushes the graphics state of the canvas on the graphics state stack. */
    void pushState() override;

    /** This function pops the graphics state of the canvas from the graphics state stack,
        restoring the most recently pushed state. */
    void popState() override;

    /** This function sets the color in the graphics state. The color is specified as three RGB
        components ranging from 0 (no intensity) to 1 (maximum intensity). If a component value is
        out of range, the result is undefined. */
    void setColor(double r, double g, double b) override;

    /** This function sets the line width in the graphics state. The line width is specified in
        units where the canvas size equals 1, and should be greater than 0 and smaller than 1. If
        the value is out of range, the result is undefined. */
    void setWidth(double w) override;

    /** This function draws a line on the canvas. The line is drawn using the current graphics
        state color and width from the point (x1,y1) to the point (x2,y2). The coordinates are
        specified in in a frame where the lower left corner of the canvas is at (0 m, 0 m) and the
        upper right corner of the canvas is at (1 m, 1 m). Segments of the line lying outside of
        the canvas are clipped away. */
    void drawLine(double x1, double y1, double x2, double y2) override;

    // ================== Data members ==================

private:
    // the (extended) canvas on which to paint; valid only during paintAndSave() execution
    class GraphicsStateCanvas;
    GraphicsStateCanvas* _canvas{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
