/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SHAPEITEM_HPP
#define SHAPEITEM_HPP

#include "ItemInfo.hpp"

////////////////////////////////////////////////////////////////////

/** The ShapeItem class is the abstract base class for all items in the 'shapes' program. It
    inherits from the Item class to offer discovery and resurrection capabilities. In addition, it
    offers facilities to handle drawing of shapes in coordination with its subclasses. The root of
    a run-time shape hierarchy must be an instance of the ShapeCanvas class, which inherits
    ShapeItem. */
class ShapeItem : public Item
{
    ITEM_ABSTRACT(ShapeItem, Item, "an abstract shape item")
    ITEM_END()

    // ================== Drawing ==================

public:
    /** This function causes the receiving shape and its children to be painted on the canvas held
        by the root of the run-time shape hierarchy. Specifically, the function calls the
        pushState() function to save the graphics state, calls its own paintSelf() function,
        forwards the paint() message to all its children (recursively) and finally calls the
        popState() function to restore the graphics state. */
    virtual void paint();

protected:
    /** This function is invoked by the paint() function to save the graphics state of the canvas
        held by the root of the run-time shape hierarchy. The implementation in this base class
        forwards the message to its parent in the run-time hierarchy until the root shape has been
        reached. */
    virtual void pushState();

    /** This function is invoked by the paint() function to restore the graphics state of the
        canvas held by the root of the run-time shape hierarchy to the most recently saved state.
        The implementation in this base class forwards the message to its parent in the run-time
        hierarchy until the root shape has been reached. */
    virtual void popState();

    /** This function must be implemented by Shape subclasses to paint the appropriate marks for
        the shape on the canvas held by the root of the run-time shape hierarchy. This can be
        accomplished by calling one or more of the drawing functions offered by this class. The
        implementation in this base class does nothing, which may be appropriate for subclasses
        such as transparent groups. */
    virtual void paintSelf();

    /** This function can be invoked from a subclass to set the color in the graphics state. The
        color is specified as three RGB components ranging from 0 (no intensity) to 1 (maximum
        intensity). If a component value is out of range, the result is undefined. The
        implementation in this base class forwards the message to its parent in the run-time
        hierarchy until the root shape has been reached. */
    virtual void setColor(double r, double g, double b);

    /** This function can be invoked from a subclass to set the line width in the graphics state.
        The line width is specified in units where the canvas size equals 1, and should be greater
        than 0 and smaller than 1. If the value is out of range, the result is undefined. The
        implementation in this base class forwards the message to its parent in the run-time
        hierarchy until the root shape has been reached. */
    virtual void setWidth(double w);

    /** This function can be invoked from a subclass to draw a line on the canvas held by the root
        of the run-time shape hierarchy. The line is drawn using the current graphics state color
        and width from the point (x1,y1) to the point (x2,y2). The coordinates are specified in in
        a frame where the lower left corner of the canvas is at (0 m, 0 m) and the upper right
        corner of the canvas is at (1 m, 1 m). Segments of the line lying outside of the canvas are
        clipped away. The implementation in this base class forwards the message to its parent in
        the run-time hierarchy until the root shape has been reached. */
    virtual void drawLine(double x1, double y1, double x2, double y2);
};

////////////////////////////////////////////////////////////////////

#endif
