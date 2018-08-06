/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTEXTINCTIONMIXINTERFACE_HPP
#define DUSTEXTINCTIONMIXINTERFACE_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** The DustExtinctionMixInterface interface offers functions to obtain the properties of a dust
    mix that are needed to calculate the effects of absorption and scattering while propagating a
    photon packet through the material. These properties include the absorption and scattering
    cross section per hydrogen atom \f$\varsigma_{\lambda}^{\text{abs}}\f$ and
    \f$\varsigma_{\lambda}^{\text{sca}}\f$, and the corresponding mass coefficients
    \f$\kappa_\lambda^{\text{abs}}\f$ and \f$\kappa_\lambda^{\text{sca}}\f$.

    It is straightforward to obtain the medium opacity from these mass coefficients and the
    medium's mass density through \f$\kappa_\lambda\rho\f$. Conversion between cross sections
    \f$\varsigma\f$ and mass coefficients \f$\kappa\f$ (through \f$\kappa=\varsigma/\mu\f$)
    requires knowledge of the mass per hydrogen atom \f$\mu\f$ of the material mix. This quantity
    is not exposed by this interface but may be used internally to perform the conversion. */
class DustExtinctionMixInterface
{
    //======== Construction - Destruction =======

protected:
    /** The empty constructor for the interface. */
    DustExtinctionMixInterface() { }

public:
    /** The empty destructor for the interface. */
    virtual ~DustExtinctionMixInterface() { }

    //======== Getters for Material Properties =======

public:
    /** This function returns the absorption cross section per hydrogen atom
        \f$\varsigma^{\text{abs}}_{\lambda}\f$ of the dust mix at wavelength \f$\lambda\f$. */
    virtual double sigmaabs(double lambda) const = 0;

    /** This function returns the scattering cross section per hydrogen atom
        \f$\varsigma^{\text{sca}}_{\lambda}\f$ of the dust mix at wavelength \f$\lambda\f$. */
    virtual double sigmasca(double lambda) const = 0;

    /** This function returns the total extinction cross section per hydrogen atom
        \f$\varsigma^{\text{ext}}_{\lambda} = \varsigma^{\text{abs}}_{\lambda} +
        \varsigma^{\text{sca}}_{\lambda}\f$ of the dust mix at wavelength \f$\lambda\f$. */
    virtual double sigmaext(double lambda) const = 0;

    /** This function returns the absorption coefficient \f$\kappa^{\text{abs}}_\lambda\f$ of the
        dust mix at wavelength \f$\lambda\f$. */
    virtual double kappaabs(double lambda) const = 0;

    /** This function returns the scattering coefficient \f$\kappa^{\text{sca}}_\lambda\f$ of the
        dust mix at wavelength \f$\lambda\f$. */
    virtual double kappasca(double lambda) const = 0;

    /** This function returns the total extinction coefficient \f$\kappa^{\text{ext}}_\lambda =
        \kappa^{\text{abs}}_\lambda + \kappa^{\text{sca}}_\lambda\f$ of the dust mix at wavelength
        \f$\lambda\f$. */
    virtual double kappaext(double lambda) const = 0;

    /** This function returns the scattering albedo \f$\varpi_\lambda =
        \varsigma_{\lambda}^{\text{sca}} / \varsigma_{\lambda}^{\text{ext}} =
        \kappa_{\lambda}^{\text{sca}} / \kappa_{\lambda}^{\text{ext}}\f$ of the dust mix at
        wavelength \f$\lambda\f$. */
    virtual double albedo(double lambda) const = 0;
};

////////////////////////////////////////////////////////////////////

#endif
