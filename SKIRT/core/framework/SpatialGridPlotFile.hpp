/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPATIALGRIDPLOTFILE_HPP
#define SPATIALGRIDPLOTFILE_HPP

#include "TextOutFile.hpp"

////////////////////////////////////////////////////////////////////

/** This class inherits from the TextOutFile class and is specifically used to write geometric
    information about a spatial grid to a data file in a text format that can be easily plotted.
    There are two format variations for 2D and 3D information, respectively. The 2D format
    describes the intersection of a spatial grid with one of the coordinate planes. The 3D format
    fully describes all or part of the spatial cells in the grid. Each line in the file contains
    two (2D) or three (3D) coordinates seperated by whitespace, or is empty. Consecutive nonempty
    lines represent a sequence of "lineto" commands; an empty line marks a "moveto" command.

    The Units object cached by the base class constructor is used by all functions to convert the
    provided coordinates to output units. */
class SpatialGridPlotFile : TextOutFile
{
    // ----------- Construction ----------

public:
    /** The constructor invokes the base class constructor to create an output file with the
        specified name and sets the appropriate precision for the numerical values in the text
        file. */
    SpatialGridPlotFile(const SimulationItem* item, string filename);

    // ----------- 2D functions ----------

public:
    /** This function outputs the specified 2D line segment. */
    void writeLine(double beg1, double beg2, double end1, double end2);

    /** This function outputs the 2D line segments describing the specified rectangle. */
    void writeRectangle(double min1, double min2, double max1, double max2);

    /** This function outputs a number of 2D line segments describing a circle with the specified
        radius centered at the specified coordinates. */
    void writeCircle(double x, double y, double radius);

    /** This function outputs a number of 2D line segments describing a circle centered at the origin with the specified
        radius. */
    void writeCircle(double radius);

    /** This function outputs a number of 2D line segments describing an arc centered at the origin with the specified
        radius and begin/end angles. */
    void writeArc(double radius, double begAngle, double endAngle);

    // ----------- 3D functions ----------

public:
    /** This function outputs the specified 3D line segment. */
    void writeLine(double x1, double y1, double z1, double x2, double y2, double z2);

    /** This function outputs 3D line segments describing the specified cuboid. */
    void writeCube(double x1, double y1, double z1, double x2, double y2, double z2);

    /** This function outputs a number of 3D line segments describing a circle with the specified
        radius at the specified height. */
    void writeCircle(double radius, double z);

    /** This function outputs a number of 3D line segments describing half a circle with the
        specified radius in the meridional plane with the specified angle. */
    void writeMeridionalHalfCircle(double radius, double phi);

    /** This function outputs a number of 3D line segments describing a coarse wireframe
        representation of a sphere with the specified center and radius, consisting of three
        orthogonal circles (one each in the planes parallel to the xy, xz, and yz coordinate planes
        through the center). The circles are drawn with fewer segments than writeCircle(), since
        this function is intended for plotting many small spheres rather than one accurate one. */
    void writeSphere(double x, double y, double z, double radius);

    /** This function outputs 3D line segments describing the specified polyhedron. Assuming
        the polyhedron has \f$n\f$ vertices, the first vector contains the \f$3n\f$ vertex
        coordinates as a sequence \f[x_0,y_0,z_0,x_1,y_1,z_1,...,x_{n-1},y_{n-1},y_{n-1}\f] The
        second vector contains the vertex indices for all polyhedron faces. Each vertex index is a
        number in the range \f$[0,n-1]\f$. The indices for each face are preceded by the number of
        vertices in the face. Assuming the polyhedron has \f$k\f$ faces, the index sequence looks
        like
        \f[m_0,i_{0,0},i_{0,1},...,i_{0,m_0-1},...,m_k,i_{k-1,0},i_{k-1,1},...,i_{k-1,m_k-1}\f] */
    void writePolyhedron(const vector<double>& coords, const vector<int>& indices);
};

////////////////////////////////////////////////////////////////////

#endif
