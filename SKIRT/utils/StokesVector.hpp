/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STOKESVECTOR_HPP
#define STOKESVECTOR_HPP

#include "Direction.hpp"

//////////////////////////////////////////////////////////////////////

/** An object of the StokesVector class describes the polarization state of a photon packet, and
    offers functions to apply certain transformations to it. Specifically, a Stokes vector contains
    the four Stokes parameters \f$I\f$, \f$Q\f$, \f$U\f$, and \f$V\f$ that are defined with respect
    to a reference direction. For efficiency reasons we store the normal to the propagation
    direction of the photon and the Stokes vector reference direction instead of the reference
    direction itself. The first parameter (\f$I\f$) indicates the total intensity, and the
    remaining parameters (\f$Q\f$, \f$U\f$, and \f$V\f$) represent various degrees of linear and
    circular polarization. Stokes parameters have a dimension of intensity. However, in our
    implementation, the parameters are normalized to dimensionless values through division by
    \f$I\f$. Consequently \f$I=1\f$ at all times, so that this parameter does not need to be
    stored. As the StokesVector does not store the propagation direction, it has to be provided
    for most transformations. The StokesVector is initialized in an unpolarized state. While
    unpolarized, the stored normal is the zero vector.*/
class StokesVector
{
public:
    // -------- constructors and setters ----------

    /** The default constructor initializes the Stokes vector to an unpolarized state. */
    StokesVector() : _polarized(false), _Q(0), _U(0), _V(0), _normal(0, 0, 0) {}

    /** This constructor initializes the Stokes vector to the specified parameter values, after
        normalizing them through division by \f$I\f$. If \f$I=0\f$, the Stokes vector is set to an
        unpolarized state. */
    StokesVector(double I, double Q, double U, double V, Direction n);

    /** This function sets the Stokes vector to an unpolarized state. */
    void setUnpolarized()
    {
        _polarized = false;
        _Q = 0;
        _U = 0;
        _V = 0;
        _normal.set(0, 0, 0);
    }

    /** This function sets the Stokes vector to the specified parameter values, after normalizing
        them through division by \f$I\f$. If \f$I=0\f$, the Stokes vector is set to an unpolarized
        state. */
    void setPolarized(double I, double Q, double U, double V, Direction n);

    /** This function copies the specified Stokes vector into the receiver. */
    void setPolarized(const StokesVector& polarization);

    // ------------ getters --------------

    /** This function returns the Stokes parameter \f$I\f$, which is always equal to one. */
    double stokesI() const { return 1.; }

    /** This function returns the Stokes parameter \f$Q\f$. */
    double stokesQ() const { return _Q; }

    /** This function returns the Stokes parameter \f$U\f$. */
    double stokesU() const { return _U; }

    /** This function returns the Stokes parameter \f$V\f$. */
    double stokesV() const { return _V; }

    /** This function returns the Stokes parameters in the provided arguments. */
    void stokes(double& I, double& Q, double& U, double& V)
    {
        I = 1.;
        Q = _Q;
        U = _U;
        V = _V;
    }

    /** This function returns the normal to the propagation direction of the photon and the
        reference direction in which the Stokes vector is defined. If the Stokes vector is in an
        unpolarized state, the zero vector is returned. */
    Direction normal() const { return _normal; }

    /** This function returns true if the Stokes vector represents polarised radiation, and false
        if it represents unpolarised radiation. */
    bool isPolarized() const { return _polarized; }

    // -------- calculated properties -------

    /** This function returns the total polarization degree for the Stokes vector. */
    double totalPolarizationDegree() const;

    /** This function returns the linear polarization degree for the Stokes vector. */
    double linearPolarizationDegree() const;

    /** This function returns the polarization position angle in radians for the Stokes vector. */
    double polarizationAngle() const;

    // ------------ transformations -----------

    /** This function adjusts the Stokes vector for a rotation of the reference axis about the
        given flight direction \f$\bf{k}\f$ over the specified angle \f$\phi\f$, clockwise when
        looking along \f$\bf{k}\f$.*/
    void rotateStokes(double phi, Direction k);

    /** This function adjusts the stokes vector reference axis so it is in the plane of the given
        propagation direction \f$\bf{k}\f$ and any direction \f$\bf{knew}\f$. The stored Stokes
        parameters are updated and the reference direction created if necessary.
        The function returns the angle by which the reference axis was rotated. */
    double rotateIntoPlane(Direction k, Direction knew);

    /** This function transforms the polarization state described by this Stokes vector by applying
        the Mueller matrix with the specified coefficients (and zero elements elsewhere) to its
        existing state. */
    void applyMueller(double S11, double S12, double S33, double S34);

private:
    bool _polarized;
    double _Q, _U, _V;
    Direction _normal;
};

//////////////////////////////////////////////////////////////////////

#endif
