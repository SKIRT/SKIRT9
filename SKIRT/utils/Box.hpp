/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BOX_HPP
#define BOX_HPP

#include "Vec.hpp"

//////////////////////////////////////////////////////////////////////

/** Box is a low-level class for working with three-dimensional "boxes": each instance represents a
    cuboid that is lined up with the cartesian coordinate axes. A box is represented by means of
    its six cartesian coordinates \f$(x_\text{min},y_\text{min},z_\text{min},
    \,x_\text{max},y_\text{max},z_\text{max})\f$ which are stored in data members with similar
    names. The class offers functions to retrieve various (mostly trivial) properties of the box,
    including its coordinates and widths in each direction and its volume. A Box instance is
    essentially immutable: once created it can no longer be changed. There is one exception to this
    rule: a derived class can replace the complete Box contents through the setExtent() function.

    The Box class is largely implemented inline (in this header file). Most compilers optimize away
    all overhead so that using this class is just as efficient as directly writing the code in
    terms of the box components. */
class Box
{
public:
    /** The default constructor creates an empty box at the origin, i.e. it initializes all box
        coordinates to zero. */
    inline Box() : _xmin(0), _ymin(0), _zmin(0), _xmax(0), _ymax(0), _zmax(0) {}

    /** This constructor initializes the box coordinates to the values provided as arguments. */
    inline Box(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax)
        : _xmin(xmin), _ymin(ymin), _zmin(zmin), _xmax(xmax), _ymax(ymax), _zmax(zmax)
    {}

    /** This constructor initializes the box coordinates to the components of the provided Vec
        objects, respectively specifying the minimum and maximum position. */
    inline Box(Vec rmin, Vec rmax)
        : _xmin(rmin.x()), _ymin(rmin.y()), _zmin(rmin.z()), _xmax(rmax.x()), _ymax(rmax.y()), _zmax(rmax.z())
    {}

    /** This function returns a reference to the receiving box object. It is useful for explicitly
        retrieving the box object from instances of classes based on Box. */
    inline const Box& extent() const { return *this; }

    /** This function stores the \f$(x_\text{min},y_\text{min},z_\text{min},
        \,x_\text{max},y_\text{max},z_\text{max})\f$ coordinates of the box in the provided
        arguments. */
    inline void extent(double& xmin, double& ymin, double& zmin, double& xmax, double& ymax, double& zmax) const
    {
        xmin = _xmin;
        ymin = _ymin;
        zmin = _zmin;
        xmax = _xmax;
        ymax = _ymax;
        zmax = _zmax;
    }

    /** This function returns the minimum position \f${\bf{r}}_\text{min}\f$ of the box. */
    inline Vec rmin() const { return Vec(_xmin, _ymin, _zmin); }

    /** This function returns the maximum position \f${\bf{r}}_\text{max}\f$ of the box. */
    inline Vec rmax() const { return Vec(_xmax, _ymax, _zmax); }

    /** This function returns the \f$x_\text{min}\f$ coordinate of the box. */
    inline double xmin() const { return _xmin; }

    /** This function returns the \f$y_\text{min}\f$ coordinate of the box. */
    inline double ymin() const { return _ymin; }

    /** This function returns the \f$z_\text{min}\f$ coordinate of the box. */
    inline double zmin() const { return _zmin; }

    /** This function returns the \f$x_\text{max}\f$ coordinate of the box. */
    inline double xmax() const { return _xmax; }

    /** This function returns the \f$y_\text{max}\f$ coordinate of the box. */
    inline double ymax() const { return _ymax; }

    /** This function returns the \f$z_\text{max}\f$ coordinate of the box. */
    inline double zmax() const { return _zmax; }

    /** This function returns the widths \f$(x_\text{max}-x_\text{min},
        y_\text{max}-y_\text{min}, z_\text{max}-z_\text{min})\f$ of the box as a Vec object. */
    inline Vec widths() const { return Vec(_xmax - _xmin, _ymax - _ymin, _zmax - _zmin); }

    /** This function returns the width \f$x_\text{max}-x_\text{min}\f$ of the box. */
    inline double xwidth() const { return _xmax - _xmin; }

    /** This function returns the width \f$y_\text{max}-y_\text{min}\f$ of the box. */
    inline double ywidth() const { return _ymax - _ymin; }

    /** This function returns the width \f$z_\text{max}-z_\text{min}\f$ of the box. */
    inline double zwidth() const { return _zmax - _zmin; }

    /** This function returns true if the position \f${\bf{r}}\f$ is inside the box, false
        otherwise. */
    inline bool contains(Vec r) const
    {
        return r.x() >= _xmin && r.x() <= _xmax && r.y() >= _ymin && r.y() <= _ymax && r.z() >= _zmin && r.z() <= _zmax;
    }

    /** This function returns true if the position \f$(x,y,z)\f$ is inside the box, false
        otherwise. */
    inline bool contains(double x, double y, double z) const
    {
        return x >= _xmin && x <= _xmax && y >= _ymin && y <= _ymax && z >= _zmin && z <= _zmax;
    }

    /** This function returns the volume \f$(x_\text{max}-x_\text{min}) \times
        (y_\text{max}-y_\text{min}) \times (z_\text{max}-z_\text{min})\f$ of the box. */
    inline double volume() const { return (_xmax - _xmin) * (_ymax - _ymin) * (_zmax - _zmin); }

    /** This function returns the length of the diagonal of the box, i.e. \f$\sqrt{
        (x_\text{max}-x_\text{min})^2 + (y_\text{max}-y_\text{min})^2 +
        (z_\text{max}-z_\text{min})^2 }\f$. */
    inline double diagonal() const
    {
        return sqrt((_xmax - _xmin) * (_xmax - _xmin) + (_ymax - _ymin) * (_ymax - _ymin)
                    + (_zmax - _zmin) * (_zmax - _zmin));
    }

    /** This function returns the position corresponding to the center of the box
        \f[ x=\tfrac12 (x_\text{min}+x_\text{max}) \f]
        \f[ y=\tfrac12 (y_\text{min}+y_\text{max}) \f]
        \f[ z=\tfrac12 (x_\text{min}+z_\text{max}) \f] */
    inline Vec center() const { return Vec(0.5 * (_xmin + _xmax), 0.5 * (_ymin + _ymax), 0.5 * (_zmin + _zmax)); }

    /** This function returns a position in the box determined by a given fraction in each spatial direction
        \f[ x=x_\text{min}+x_\text{frac}\times(x_\text{max}-x_\text{min}) \f]
        \f[ y=y_\text{min}+y_\text{frac}\times(y_\text{max}-y_\text{min}) \f]
        \f[ z=z_\text{min}+z_\text{frac}\times(z_\text{max}-z_\text{min}) \f]
        The specified fractions must be between zero and one; this is \em not checked by the function. */
    inline Vec fracPos(double xfrac, double yfrac, double zfrac) const
    {
        return Vec(_xmin + xfrac * (_xmax - _xmin), _ymin + yfrac * (_ymax - _ymin), _zmin + zfrac * (_zmax - _zmin));
    }

    /** This function returns a position in the box determined by a given fraction in each spatial direction.
        \f[ x=x_\text{min}+\frac{x_\text{d}}{x_\text{n}}\times(x_\text{max}-x_\text{min}) \f]
        \f[ y=y_\text{min}+\frac{y_\text{d}}{y_\text{n}}\times(y_\text{max}-y_\text{min}) \f]
        \f[ z=z_\text{min}+\frac{z_\text{d}}{z_\text{n}}\times(z_\text{max}-z_\text{min}) \f]
        Each of the fractions is specified as the quotient of two integers; the integers are
        converted to floating point before the division is made. The quotients must be between zero
        and one; this is \em not checked by the function. */
    inline Vec fracPos(int xd, int yd, int zd, int xn, int yn, int zn) const
    {
        return Vec(_xmin + xd * (_xmax - _xmin) / xn, _ymin + yd * (_ymax - _ymin) / yn,
                   _zmin + zd * (_zmax - _zmin) / zn);
    }

    /** This function calculates the cell indices for a given position, assuming that the box would
        be partitioned in a given number of cells in each spatial direction. */
    inline void cellIndices(int& i, int& j, int& k, Vec r, int nx, int ny, int nz) const
    {
        i = std::max(0, std::min(nx - 1, static_cast<int>(nx * (r.x() - _xmin) / (_xmax - _xmin))));
        j = std::max(0, std::min(ny - 1, static_cast<int>(ny * (r.y() - _ymin) / (_ymax - _ymin))));
        k = std::max(0, std::min(nz - 1, static_cast<int>(nz * (r.z() - _zmin) / (_zmax - _zmin))));
    }

    /** This function intersects the receiving axis-aligned bounding box with a ray (half-line)
        defined by the specified starting position \f$\bf{r}\f$ and direction \f$\bf{k}\f$. If the
        ray intersects the box, the function returns true after storing the near and far
        intersection distances relative to the starting position in \f$s_\mathrm{min}\f$ and
        \f$s_\mathrm{max}\f$, respectively. If the starting position is inside the box, the nearest
        intersection distance is set to zero. In other words, the following relation always holds:
        \f$0\le s_\mathrm{min} < s_\mathrm{max}\f$.

        If the ray does not intersect the box, the function returns false and the values of
        \f$s_\mathrm{min}\f$ and \f$s_\mathrm{max}\f$ are undefined. A ray that "touches" the box
        border in a single point is not considered to intersect the box. A ray along an edge or
        face on the lower side of the box is considered to intersect the box, while ray along an
        edge or face on the upper side of the box is considered \em not to intersect the box. This
        approach avoids duplicate intersection of adjacent boxes.

        The function employs the slab method originated by Kay and Kajiya (1986) and adapted by
        Haines (1989) as described in "Geometric Tools for Computer Graphics" by Scheider and
        Eberly (2003, Elsevier). */
    bool intersects(Vec r, const Vec k, double& smin, double& smax) const;

    /** This function intersects the receiving axis-aligned bounding box with a sphere defined by
        the specified center position \f${\bf{r}}_\mathrm{c}\f$ and radius \f$r\f$. It returns true
        if the box and the sphere intersect, and false otherwise.

        The function employs the algorithm due to Jim Arvo described in "Graphics Gems" (1990). */
    bool intersects(Vec rc, double r) const;

protected:
    /** This function replaces the extent of the box with the newly specified values. This function
        is intended for use in derived classes only. */
    void setExtent(const Box& extent) { *this = extent; }

    /** This function replaces the extent of the box with the newly specified values. This function
        is intended for use in derived classes only. */
    void setExtent(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax)
    {
        _xmin = xmin;
        _ymin = ymin;
        _zmin = zmin;
        _xmax = xmax;
        _ymax = ymax;
        _zmax = zmax;
    }

private:
    /** These data members represent the cartesian vector components */
    double _xmin, _ymin, _zmin, _xmax, _ymax, _zmax;
};

//////////////////////////////////////////////////////////////////////

#endif
