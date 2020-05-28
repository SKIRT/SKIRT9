/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef COSMOLOGY_HPP
#define COSMOLOGY_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** An instance of a Cosmology subclass specifies the cosmology related parameters of the simulated
    model, including the redshift at which the model coordinate frame is assumed to reside, as well
    as the spatial geometry and expansion behavior of the universe. */
class Cosmology : public SimulationItem
{
    ITEM_ABSTRACT(Cosmology, SimulationItem, "a set of cosmology parameters, including redshift")
    ITEM_END()

    //============== Functions to be implemented in subclasses ============

public:
    /** This function returns the redshift at which the model resides, or zero if it resides in the
        Local Universe. */
    virtual double modelRedshift() const = 0;

    /** This function returns the angular-diameter distance \f$d_\mathrm{A}(z)\f$, which converts a
        \em proper transverse separation \f$\mathrm{d}l\f$ to the corresponding observed angular
        separation \f$\mathrm{d}\psi\f$, \f[ \mathrm{d}\psi = \frac{\mathrm{d}l}{d_\mathrm{A}(z)}.
        \f] */
    virtual double angularDiameterDistance() const = 0;

    /** This function returns the luminosity distance, which converts a total luminosity
        \f$L_\mathrm{tot}\f$ to the corresponding observed total flux \f$F_\mathrm{tot}\f$, \f[
        F_\mathrm{tot} = \frac{L_\mathrm{tot}}{4\pi\,d_\mathrm{L}^2(z)} \f] or a neutral-style
        monochromatic luminosity \f$\lambda L_\lambda\f$ to the corresponding observed
        neutral-style flux, \f[ (1+z)\lambda \, F_\lambda[(1+z)\lambda] = \frac{\lambda
        L_\lambda[\lambda]}{4\pi\,d_\mathrm{L}^2(z)} \f] where \f$\lambda\f$ is the emitted
        wavelength and \f$(1+z)\lambda\f$ is the observed wavelength. */
    virtual double luminosityDistance() const = 0;

    /** This function returns the relative expansion rate of the universe, \f$\dot{a}/a\f$, where
        \f$a\f$ is the dimensionless scale factor. The relative expansion rate has units of
        \f$\mathrm{s}^{-1}\f$. */
    virtual double relativeExpansionRate() const = 0;
};

////////////////////////////////////////////////////////////////////

#endif
