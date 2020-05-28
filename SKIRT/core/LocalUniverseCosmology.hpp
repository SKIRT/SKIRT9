/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LOCALUNIVERSECOSMOLOGY_HPP
#define LOCALUNIVERSECOSMOLOGY_HPP

#include "Cosmology.hpp"

////////////////////////////////////////////////////////////////////

/** An instance of the LocalUniverseCosmology class represents a special-case "cosmology" where the
    simulated model resides at redshift zero. As a result, there is no need to specify the spatial
    geometry and expansion behavior of the universe. */
class LocalUniverseCosmology : public Cosmology
{
    ITEM_CONCRETE(LocalUniverseCosmology, Cosmology, "the model is at redshift zero in the Local Universe")
    ITEM_END()

    //============== Functions overridden from base class ============

public:
    /** This function returns a zero redshift, indicating that the model resides in the Local
        Universe. */
    double modelRedshift() const override;

    /** This function returns the angular-diameter distance \f$d_\mathrm{A}(z)\f$. For this class,
        it always returns zero, indicating that the model resides in the Local Universe. */
    double angularDiameterDistance() const override;

    /** This function returns the luminosity distance. For this class, it always returns zero,
        indicating that the model resides in the Local Universe. */
    double luminosityDistance() const override;

    /** This function returns the present-day relative expansion rate of the universe, often called
        the Hubble constant, assuming a fixed value of 67.5 (km/s)/Mpc. The function returns this
        value in SI units, i.e. \f$\mathrm{s}^{-1}\f$.

        To configure another value of the present-day relative expansion rate, use the
        FlatUniverseCosmology class with a redshift of zero and the desired value for the reduced
        Hubble parameter. */
    double relativeExpansionRate() const override;
};

////////////////////////////////////////////////////////////////////

#endif
