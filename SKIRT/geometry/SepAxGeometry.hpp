/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SEPAXGEOMETRY_HPP
#define SEPAXGEOMETRY_HPP

#include "AxGeometry.hpp"

//////////////////////////////////////////////////////////////////////

/** The SepAxGeometry class is an abstract subclass of the AxGeometry class, and
    serves as a base class for axisymmetric geometries, where the density is a separable
    function of \f$R\f$ and \f$z\f$. This means that we can write it as \f[ \rho({\bf{r}}) =
    \rho_R(R)\,\rho_z(z), \f] with \f$\rho_R\f$ and \f$\rho_z\f$ functions that satisfy the
    normalization \f[ 2\pi \int_0^\infty \rho_R(R)\, R\, {\text{d}}R = \int_{-\infty}^\infty
    \rho_z(z)\, {\text{d}}z = 1. \f] */
class SepAxGeometry : public AxGeometry
{
    ITEM_ABSTRACT(SepAxGeometry, AxGeometry, "a separable axisymmetric geometry")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function generates a random position from the geometry, by drawing a random
        point from the three-dimensional probability density \f$p({\bf{r}})\, {\text{d}}{\bf{r}} =
        \rho({\bf{r}})\, {\text{d}}{\bf{r}}\f$. In a separable axisymmetric geometry, this
        three-dimensional probability density is separable into a cylindrical, an azimuthal and a
        vertical component: \f[ p({\bf{r}})\, {\text{d}}{\bf{r}} = \left[
        2\pi\,\rho_R(R)\,R\,{\text{d}}R \right] \left[ \frac{{\text{d}}\phi}{2\pi} \right] \left[
        \rho_z(z)\,{\text{d}}z \right]. \f] A random position can hence be constructed by combining
        random cylindrical coordinates, each chosen from their own probability distributions. A
        random azimuth \f$\phi\f$ is readily found by chosing a random deviate \f${\cal{X}}\f$ and
        setting \f$ \phi = 2\pi {\cal{X}} \f$. A random cylindrical radius and height can be found
        through the member functions randomR() and randomz(), which have to be provided for each
        class derived from the SepAxGeometry class. */
    Position generatePosition() const override;

    /** This pure virtual function returns the cylindrical radius \f$R\f$ of a random position. For
        any axisymmetric geometry in which the density is separable in \f$R\f$ and \f$z\f$, this
        can be achieved by generating a random \f$R\f$ from the one-dimensional probability
        distribution \f[ p(R)\, {\text{d}}R = 2\pi\, \rho_R(R)\, R\, {\text{d}}R. \f] */
    virtual double randomCylRadius() const = 0;

    /** This pure virtual function returns the height \f$z\f$ of a random position. For any
        axisymmetric geometry in which the density is separable in \f$R\f$ and \f$z\f$, this can be
        achieved by generating a random \f$z\f$ from the one-dimensional probability distribution
        \f[ p(z)\, {\text{d}}z = \rho_z(z)\, {\text{d}}z. \f] */
    virtual double randomZ() const = 0;
};

////////////////////////////////////////////////////////////////////

#endif
