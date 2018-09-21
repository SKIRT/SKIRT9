/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef AXGEOMETRY_HPP
#define AXGEOMETRY_HPP

#include "Geometry.hpp"

//////////////////////////////////////////////////////////////////////

/** The AxGeometry class is an abstract subclass of the general Geometry class, and
    serves as a base class for axisymmetric geometries, i.e. geometries where the
    density can be written as \f$\rho({\bf{r}}) = \rho(R,z)\f$. */
class AxGeometry : public Geometry
{
    ITEM_ABSTRACT(AxGeometry, Geometry, "an axisymmetric geometry")
        ATTRIBUTE_TYPE_INSERT(AxGeometry, "Dimension2")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the geometry, which is 2 for all subclasses
        of this class since it is a base class for axisymmetric geometries. */
    int dimension() const override;

    /** This function returns the density \f$\rho({\bf{r}})\f$ at the position \f${\bf{r}}\f$. For
        axisymmetric geometries, the density only depends on the cylindrical coordinates \f$R\f$
        and \f$z\f$, such that only the member function \f$\rho(R,z)\f$ needs to be defined. */
    double density(Position bfr) const override;

    /** This pure virtual function returns the density \f$\rho(R,z)\f$ at the cylindrical radius
        \f$R\f$ and height \f$z\f$. */
    virtual double density(double R, double z) const = 0;

    /** This function returns the X-axis surface density, i.e. the integration of the density
        along the entire X-axis, \f[ \Sigma_X = \int_{-\infty}^\infty \rho(x,0,0)\,{\text{d}}x. \f]
        For an axisymmetric geometry, we can write \f$ \Sigma_X = 2\,\Sigma_R \f$, where \f$
        \Sigma_R \f$ is the radial surface density. */
    double SigmaX() const override;

    /** This function returns the Y-axis surface density, i.e. the integration of the density
        along the entire Y-axis, \f[ \Sigma_Y = \int_{-\infty}^\infty \rho(0,y,0)\,{\text{d}}y. \f]
        For an axisymmetric geometry, we can write \f$ \Sigma_Y = 2\,\Sigma_R \f$, where \f$
        \Sigma_R \f$ is the radial surface density. */
    double SigmaY() const override;

    /** This pure virtual function returns the radial surface density, i.e. the integration of
        the density along a line in the equatorial plane starting at the centre of the coordinate
        system, \f[ \Sigma_R = \int_0^\infty \rho(R,0)\,{\text{d}}R. \f] */
    virtual double SigmaR() const = 0;
};

////////////////////////////////////////////////////////////////////

#endif
