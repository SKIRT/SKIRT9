/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CYLCELL_HPP
#define CYLCELL_HPP

#include "Box.hpp"

//////////////////////////////////////////////////////////////////////

/** CylCell is a low-level class for working with basic three-dimensional cells in cylindrical
    coordinates. Each CylCell instance represents a cell bordered by:

    - two concentric vertical cylinders defined by \f$ 0 \le R_\text{min} \le R_\text{max} \f$,

    - two meridional half-planes (with \f$R>0\f$) defined by \f$ -\pi \le \varphi_\text{min} \le
    \varphi_\text{max} \le \pi \f$,

    - two horizontal planes defined by \f$ z_\text{min} \le z_\text{max}\f$.

    These six values are stored in data members with similar names during construction.

    \note Because of the limitations on the range of \f$\varphi\f$, it is not possible for a
    Cylcell to straddle the negative x-axis of the Cartesian model coordinate system.

    The class offers functions to retrieve various basic properties of the cell, such as its border
    coordinates and its volume, and for geometric operations such as determining whether a given
    Cartesian position is inside the cell. */
class CylCell
{
public:
    /** The default constructor creates an empty cell at the origin, i.e. it initializes all border
        coordinates to zero. */
    CylCell() : _Rmin(0), _phimin(0), _zmin(0), _Rmax(0), _phimax(0), _zmax(0) {}

    /** This constructor initializes the cell border coordinates to the values provided as
        arguments. It does not verify that these values conform to the limits described in the
        class header. Non-comforming values lead to undefined behavior. */
    CylCell(double Rmin, double phimin, double zmin, double Rmax, double phimax, double zmax)
        : _Rmin(Rmin), _phimin(phimin), _zmin(zmin), _Rmax(Rmax), _phimax(phimax), _zmax(zmax)
    {}

    /** This function stores the cell border coordinates in the provided arguments. */
    void extent(double& Rmin, double& phimin, double& zmin, double& Rmax, double& phimax, double& zmax) const
    {
        Rmin = _Rmin;
        phimin = _phimin;
        zmin = _zmin;
        Rmax = _Rmax;
        phimax = _phimax;
        zmax = _zmax;
    }

    /** This function returns the \f$R_\text{min}\f$ border coordinate of the cell. */
    double Rmin() const { return _Rmin; }

    /** This function returns the \f$\varphi_\text{min}\f$ border coordinate of the cell. */
    double phimin() const { return _phimin; }

    /** This function returns the \f$z_\text{min}\f$ border coordinate of the cell. */
    double zmin() const { return _zmin; }

    /** This function returns the \f$R_\text{max}\f$ border coordinate of the cell. */
    double Rmax() const { return _Rmax; }

    /** This function returns the \f$\varphi_\text{max}\f$ border coordinate of the cell. */
    double phimax() const { return _phimax; }

    /** This function returns the \f$z_\text{max}\f$ border coordinate of the cell. */
    double zmax() const { return _zmax; }

    /** This function returns the width \f$R_\text{max}-R_\text{min}\f$ of the cell. */
    double Rwidth() const { return _Rmax - _Rmin; }

    /** This function returns the width \f$\varphi_\text{max}-\varphi_\text{min}\f$ of the cell. */
    double phiwidth() const { return _phimax - _phimin; }

    /** This function returns the width \f$z_\text{max}-z_\text{min}\f$ of the cell. */
    double zwidth() const { return _zmax - _zmin; }

    /** This function returns the volume of the cell, given by \f$\frac{1}{2}
        (R_\text{max}^2-R_\text{min}^2) (\varphi_\text{max}-\varphi_\text{min})
        (z_\text{max}-z_\text{min})\f$. */
    double volume() const { return 0.5 * (_Rmax * _Rmax - _Rmin * _Rmin) * (_phimax - _phimin) * (_zmax - _zmin); }

    /** This function returns true if the position \f$(R,\varphi,z)\f$ in cylindrical coordinates is
        inside the cell, and false otherwise. */
    bool contains(double R, double phi, double z) const
    {
        return R >= _Rmin && R <= _Rmax && phi >= _phimin && phi <= _phimax && z >= _zmin && z <= _zmax;
    }

    /** This function returns true if the Cartesian position \f${\bf{r}}=(x,y,z)\f$ is inside the
        cell, and false otherwise. */
    bool contains(Vec r) const;

    /** This function returns the Cartesian bounding box of the cell, in other words the smallest
        cuboid lined up with the Cartesian coordinate axes that encloses the cell. */
    Box boundingBox() const;

private:
    // These data members represent the cylindrical border coordinates
    double _Rmin, _phimin, _zmin, _Rmax, _phimax, _zmax;
};

//////////////////////////////////////////////////////////////////////

#endif
