/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CANVAS_HPP
#define CANVAS_HPP

#include "Table.hpp"
class ConvexPolygon;

////////////////////////////////////////////////////////////////////

/** The Canvas class represents a simple pixelated canvas with support for drawing RGB-colored
    lines and saving to a TIFF file for viewing. The canvas is always square with a logical
    coordinate range between 0 and 1 in each direction. The number of pixels must be specified at
    construction and is the same in each direction. Colors are specified as three RGB components
    ranging from 0 (no intensity) to 1 (maximum intensity). Line widths are specified in units
    where the canvas size equals 1. The canvas is initialized with a black background (i.e. all
    pixels values are set to zero), and with a default line color of white (i.e. the three
    components equal 1) and a default line width of 0.01 (i.e. 1/100 of the canvas size). */
class Canvas
{
    // ================== Constructing ==================

public:
    /** Constructs a canvas with the specified number of pixels in each direction, which must be at
        least 3 and at most 32000 (the number must fit in 15 bits). If the number of pixels is out
        of range, the behavior is undefined. The canvas is initialized with a black background,
        with a default line color of white, and a default line width of 0.01. These defaults can be
        changed with the respective setters. */
    Canvas(size_t numPixels);

    // ================== Graphics state and drawing ==================

public:
    /** Sets the color of lines drawn by subsequent invocations of drawLine(). The color is
        specified as three RGB components ranging from 0 (no intensity) to 1 (maximum intensity).
        If a component value is out of range, the result is undefined. */
    void setColor(double r, double g, double b);

    /** Sets the width of lines drawn by subsequent invocations of drawLine(). The line width is
        specified in units where the canvas size equals 1, and should be greater than 0 and smaller
        than 1. If the value is out of range, the result is undefined. */
    void setWidth(double w);

    /** Fills a rectangle lined up with the coordinate axes with corner points (x1,y1) and (x2,y2)
        using the current color. The coordinates are specified in the canvas coordinate system
        ranging from 0 and 1 in each direction. Segments of the rectangle lying outside of the
        canvas are clipped away. */
    void fillBox(double x1, double y1, double x2, double y2);

    /** Fills the specified convex polygon using the current color. The coordinates are specified
        in the canvas coordinate system ranging from 0 and 1 in each direction. Segments of the
        polygon lying outside of the canvas are clipped away. If the points of the polygon
        violate the restrictions described for the ConvexPolygon class, the behavior of this
        function is undefined. */
    void fillConvexPolygon(const ConvexPolygon& poly);

    /** Draws a line with the current color and width from the point (x1,y1) to the point (x2,y2).
        The coordinates are specified in the canvas coordinate system ranging from 0 and 1 in each
        direction. Segments of the line lying outside of the canvas are clipped away. More
        specifically, the function proceeds as follows. It first constructs a virtual rectangle by
        expanding the mathematical line between the given points with half of the current width in
        all directions. It then assigns the current color to all pixels for which the pixel center
        lies inside that rectangle (and leaves any other pixels alone). */
    void drawLine(double x1, double y1, double x2, double y2);

    // ================== Saving ==================

public:
    /** Saves the current contents of the canvas as a TIFF file with the specified file path. The
        file path should already include the appropriate filename extension. If the operation was
        successful, the function returns true; otherwise it returns false. */
    bool saveToTiff(string filepath) const;

    // ================== Private utilities ==================

private:
    /** Determines the indices for the range of pixels that lie inside the specified coordinate
        range. If the range contains at least one pixel, the function places the indices in the
        output arguments and returns true. If not, the function returns false and the contents of
        the output arguments is undefined. */
    bool indicesInRange(double x1, double x2, size_t& i1, size_t& i2);

    // ================== Data members ==================

private:
    Table<3> _pixels;  // the pixel values (3 RGB components for each pixel)
    size_t _n = 0;     // the number of pixels
    double _r = 1.;    // current color
    double _g = 1.;
    double _b = 1.;
    double _w2 = 0.005;  // half of the current width
};

////////////////////////////////////////////////////////////////////

#endif
