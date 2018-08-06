/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ELECTRONEXTINCTIONMIXINTERFACE_HPP
#define ELECTRONEXTINCTIONMIXINTERFACE_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** The ElectronExtinctionMixInterface interface offers functions to obtain the properties of an
    electron population that are needed to calculate the effects of absorption and scattering while
    propagating a photon packet through the material. These properties include the absorption and
    scattering cross section per electron \f$\varsigma_{\lambda}^{\text{abs}}\f$ and
    \f$\varsigma_{\lambda}^{\text{sca}}\f$. The absorption cross section for an electron is
    trivially zero, and the scattering cross section is wavelength independet, but the functions
    are provided in this interface by analogy to the corresponding interface for gas.

    It is straightforward to obtain the medium opacity from these cross sections and the electron
    number density \f$n\f$ through \f$\varsigma_{\lambda}n\f$. */
class ElectronExtinctionMixInterface
{
    //======== Construction - Destruction =======

protected:
    /** The empty constructor for the interface. */
    ElectronExtinctionMixInterface() { }

public:
    /** The empty destructor for the interface. */
    virtual ~ElectronExtinctionMixInterface() { }

    //======== Getters for Material Properties =======

public:
    /** This function returns the absorption cross section per electron
        \f$\varsigma^{\text{abs}}_{\lambda}\f$ of the material mix at wavelength \f$\lambda\f$. */
    virtual double sigmaabs(double lambda) const = 0;

    /** This function returns the scattering cross section per electron
        \f$\varsigma^{\text{sca}}_{\lambda}\f$ of the material mix at wavelength \f$\lambda\f$. */
    virtual double sigmasca(double lambda) const = 0;

    /** This function returns the total extinction cross section per electron
        \f$\varsigma^{\text{ext}}_{\lambda} = \varsigma^{\text{abs}}_{\lambda} +
        \varsigma^{\text{sca}}_{\lambda}\f$ of the material mix at wavelength \f$\lambda\f$. */
    virtual double sigmaext(double lambda) const = 0;

    /** This function returns the scattering albedo \f$\varpi_\lambda =
        \varsigma_{\lambda}^{\text{sca}} / \varsigma_{\lambda}^{\text{ext}} =
        \kappa_{\lambda}^{\text{sca}} / \kappa_{\lambda}^{\text{ext}}\f$ of the material mix at
        wavelength \f$\lambda\f$. */
    virtual double albedo(double lambda) const = 0;
};

////////////////////////////////////////////////////////////////////

#endif
