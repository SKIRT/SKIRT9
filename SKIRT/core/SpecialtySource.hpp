/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPECIALTYSOURCE_HPP
#define SPECIALTYSOURCE_HPP

#include "NormalizedSource.hpp"
#include "VelocityInterface.hpp"

//////////////////////////////////////////////////////////////////////

/** SpecialtySource is an abstract class representing a primary radiation source characterized by a
    single SED object and a single bulk velocity, i.e. both the spectrum and the bulk velocity are
    identical in all locations. This combination is suitable for specialty sources for which a
    velocity field makes no sense.

    Subclasses must handle the spatial distribution of the source, and can optionally add
    anisotropy and/or polarization. */
class SpecialtySource : public NormalizedSource, public VelocityInterface
{
    ITEM_ABSTRACT(SpecialtySource, NormalizedSource, "a primary source with a single bulk velocity")

        ATTRIBUTE_SUB_PROPERTIES_HERE(SpecialtySource)

        PROPERTY_DOUBLE(velocityX, "the bulk velocity of the source, x component")
        ATTRIBUTE_QUANTITY(velocityX, "velocity")
        ATTRIBUTE_MIN_VALUE(velocityX, "[-100000 km/s")
        ATTRIBUTE_MAX_VALUE(velocityX, "100000 km/s]")
        ATTRIBUTE_DEFAULT_VALUE(velocityX, "0")
        ATTRIBUTE_RELEVANT_IF(velocityX, "Panchromatic")
        ATTRIBUTE_DISPLAYED_IF(velocityX, "Level2")
        ATTRIBUTE_INSERT(velocityX, "Panchromatic&velocityX:Dimension3")

        PROPERTY_DOUBLE(velocityY, "the bulk velocity of the source, y component")
        ATTRIBUTE_QUANTITY(velocityY, "velocity")
        ATTRIBUTE_MIN_VALUE(velocityY, "[-100000 km/s")
        ATTRIBUTE_MAX_VALUE(velocityY, "100000 km/s]")
        ATTRIBUTE_DEFAULT_VALUE(velocityY, "0")
        ATTRIBUTE_RELEVANT_IF(velocityY, "Panchromatic")
        ATTRIBUTE_DISPLAYED_IF(velocityY, "Level2")
        ATTRIBUTE_INSERT(velocityY, "Panchromatic&velocityY:Dimension3")

        PROPERTY_DOUBLE(velocityZ, "the bulk velocity of the source, z component")
        ATTRIBUTE_QUANTITY(velocityZ, "velocity")
        ATTRIBUTE_MIN_VALUE(velocityZ, "[-100000 km/s")
        ATTRIBUTE_MAX_VALUE(velocityZ, "100000 km/s]")
        ATTRIBUTE_DEFAULT_VALUE(velocityZ, "0")
        ATTRIBUTE_RELEVANT_IF(velocityZ, "Panchromatic")
        ATTRIBUTE_DISPLAYED_IF(velocityZ, "Level2")
        ATTRIBUTE_INSERT(velocityZ, "Panchromatic&velocityZ:Dimension2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function creates a private object offering the velocity interface if the bulk velocity
        is nonzero. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the source, which is the same as the dimension of
        its spatial distribution (to be provided by the subclass), except if there is a nonzero
        bulk velocity. */
    int dimension() const override;

    /** This function returns true if the velocity of the source is nonzero. */
    bool hasVelocity() const override;

    /** This function implements the VelocityInterface interface. It returns the bulk velocity
         of this source, as configured by the user. */
    Vec velocity() const override;

    /** This function causes the photon packet \em pp to be launched from the source. It passes the
        request through to the subclass, adding the redshift interface for the configured bulk
        velocity to the list of arguments. */
    void launchNormalized(PhotonPacket* pp, size_t historyIndex, double lambda, double Lw) const override;

    //============== Functions to be implemented in each subclass =============

    /** This function returns the dimension of the spatial distribution implemented by the
        subclass, taking into account anisotropic emission or polarization, if any. It must be
        implemented in a subclass. */
    virtual int geometryDimension() const = 0;

    /** This function causes the photon packet \em pp to be launched from the source using the
        given history index, wavelength, weighted luminosity contribution, and redshift interface
        (corresponding to the bulk velocity of the source). It must be implemented in a subclass to
        handle the spatial distribution of the source, optionally adding anisotropy and/or
        polarization. */
    virtual void launchSpecialty(PhotonPacket* pp, size_t historyIndex, double lambda, double Lw,
                                 VelocityInterface* bvi) const = 0;

    //======================== Data Members ========================

private:
    // pointer to an object offering the redshift interface; either "this" or null pointer if the bulk velocity is zero
    VelocityInterface* _bvi{nullptr};
};

//////////////////////////////////////////////////////////////////////

#endif
