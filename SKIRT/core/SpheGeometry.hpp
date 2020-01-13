/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHEGEOMETRY_HPP
#define SPHEGEOMETRY_HPP

#include "Geometry.hpp"

////////////////////////////////////////////////////////////////////

/** The SpheGeometry class is an abstract subclass of the general Geometry class, and
    serves as a base class for spherically symmetric geometries, i.e. geometries where the
    density depends only on the modulus \f$r = |{\bf{r}}|\f$ of the position vector. */
class SpheGeometry : public Geometry
{
    ITEM_ABSTRACT(SpheGeometry, Geometry, "a spherically symmetric geometry")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the geometry, which is 1 for all subclasses of this
        class since it is a base class for spherically symmetric geometries. */
    int dimension() const override;

    /** This function returns the density \f$\rho({\bf{r}})\f$ at the position \f${\bf{r}}\f$. For
        spherically symmetric geometries, the density only depends on the modulus
        \f$r=|{\bf{r}}|\f$ of the position vector, such that only the member function \f$\rho(r)\f$
        needs to be defined. */
    double density(Position bfr) const override;

    /** This pure virtual function returns the density \f$\rho(r)\f$ at a radius \f$r\f$. */
    virtual double density(double r) const = 0;

    /** This pure virtual function generates a random position from the geometry, by
        drawing a random point from the three-dimensional probability density \f$p({\bf{r}})\,
        {\text{d}}{\bf{r}} = \rho({\bf{r}})\, {\text{d}}{\bf{r}}\f$. In spherical symmetry, this
        three-dimensional probability density is separable into a radial component and a direction
        component: \f[ p({\bf{r}})\, {\text{d}}{\bf{r}} = \left[ 4\pi\, \rho(r)\, r^2\, {\text{d}}r
        \right] \left[ \frac{{\text{d}}{\bf{k}}}{4\pi} \right]. \f] A random position can hence be
        constructed by combining a random direction on the unit sphere with a radius randomly
        generated from the one-dimensional probability density \f[ p(r)\,{\text{d}}r = 4\pi\,
        \rho(r)\, r^2\, {\text{d}}r. \f] The former is found through the general function
        randomDirection(), whereas the latter can be found through the member function
        randomradius(), which has to be provided in each class derived from the SpheGeometry
        class. */
    Position generatePosition() const override;

    /** This pure virtual function returns a radius randomly generated from the one-dimensional
        probability distribution \f[ p(r)\, {\text{d}}r = 4\pi\, \rho(r)\, r^2\, {\text{d}}r. \f]
        */
    virtual double randomRadius() const = 0;

    /** This function returns the X-axis surface density, i.e. the integration of the density
        along the entire X-axis, \f[ \Sigma_X = \int_{-\infty}^\infty \rho(x,0,0)\,{\text{d}}x. \f]
        For a spherical geometry, we can write \f$ \Sigma_X = 2\,\Sigma_r \f$, where \f$
        \Sigma_r \f$ is the radial surface density. */
    double SigmaX() const override;

    /** This function returns the Y-axis surface density, i.e. the integration of the density
        along the entire Y-axis, \f[ \Sigma_Y = \int_{-\infty}^\infty \rho(0,y,0)\,{\text{d}}y. \f]
        For a spherical geometry, we can write \f$ \Sigma_Y = 2\,\Sigma_r \f$, where \f$
        \Sigma_r \f$ is the radial surface density. */
    double SigmaY() const override;

    /** This function returns the Z-axis surface density, i.e. the integration of the density
        along the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\,{\text{d}}z. \f]
        For a spherical geometry, we can write \f$ \Sigma_Z = 2\,\Sigma_r \f$, where \f$
        \Sigma_r \f$ is the radial surface density. */
    double SigmaZ() const override;

    /** This pure virtual function returns the radial surface density, i.e. the integration of
        the density along a line starting at the centre of the coordinate
        system, \f[ \Sigma_r = \int_0^\infty \rho(r)\,{\text{d}}r. \f] */
    virtual double Sigmar() const = 0;
};

////////////////////////////////////////////////////////////////////

#endif
