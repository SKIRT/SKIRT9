/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NORMALIZEDSOURCE_HPP
#define NORMALIZEDSOURCE_HPP

#include "Source.hpp"
#include "BulkVelocityInterface.hpp"
#include "LuminosityNormalization.hpp"
#include "SED.hpp"

//////////////////////////////////////////////////////////////////////

/** NormalizedSource is an abstract class representing a primary radiation source characterized by
    a single SED object, i.e. the spectral distribution is identical in all spatial locations. The
    source can have a single bulk velocity, i.e. the bulk velocity is also identical in all
    locations. The bolometric power of the source is characterized by a LuminosityNormalization
    object.

    Subclasses must handle the spatial distribution of the source, and can optionally add
    anisotropy and/or polarization. */
class NormalizedSource : public Source, public BulkVelocityInterface
{
    ITEM_ABSTRACT(NormalizedSource, Source, "a primary source with a single SED")

    ATTRIBUTE_SUB_PROPERTIES_HERE(NormalizedSource)

    PROPERTY_ITEM(sed, SED, "the spectral energy distribution for the source")
        ATTRIBUTE_DEFAULT_VALUE(sed, "BlackBodySED")

    PROPERTY_ITEM(normalization, LuminosityNormalization, "the type of luminosity normalization for the source")
        ATTRIBUTE_DEFAULT_VALUE(normalization, "IntegratedLuminosityNormalization")

    PROPERTY_DOUBLE(velocityX, "the bulk velocity of the source, x component")
        ATTRIBUTE_QUANTITY(velocityX, "velocity")
        ATTRIBUTE_MIN_VALUE(velocityX, "[0")
        ATTRIBUTE_MAX_VALUE(velocityX, "100000 km/s]")
        ATTRIBUTE_DEFAULT_VALUE(velocityX, "0")
        ATTRIBUTE_RELEVANT_IF(velocityX, "Panchromatic")
        ATTRIBUTE_DISPLAYED_IF(velocityX, "Level2")
        ATTRIBUTE_INSERT(velocityX, "Panchromatic&velocityX:Dimension3")

    PROPERTY_DOUBLE(velocityY, "the bulk velocity of the source, y component")
        ATTRIBUTE_QUANTITY(velocityY, "velocity")
        ATTRIBUTE_MIN_VALUE(velocityY, "[0")
        ATTRIBUTE_MAX_VALUE(velocityY, "100000 km/s]")
        ATTRIBUTE_DEFAULT_VALUE(velocityY, "0")
        ATTRIBUTE_RELEVANT_IF(velocityY, "Panchromatic")
        ATTRIBUTE_DISPLAYED_IF(velocityY, "Level2")
        ATTRIBUTE_INSERT(velocityY, "Panchromatic&velocityY:Dimension3")

    PROPERTY_DOUBLE(velocityZ, "the bulk velocity of the source, z component")
        ATTRIBUTE_QUANTITY(velocityZ, "velocity")
        ATTRIBUTE_MIN_VALUE(velocityZ, "[0")
        ATTRIBUTE_MAX_VALUE(velocityZ, "100000 km/s]")
        ATTRIBUTE_DEFAULT_VALUE(velocityZ, "0")
        ATTRIBUTE_RELEVANT_IF(velocityZ, "Panchromatic")
        ATTRIBUTE_DISPLAYED_IF(velocityZ, "Level2")
        ATTRIBUTE_INSERT(velocityZ, "Panchromatic&velocityZ:Dimension2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function creates a private object offering the redshift interface if the bulk velocity
        is nonzero. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the source, which is the same as the dimension of
        its spatial distribution (to be provided by the subclass), except if there is a nonzero
        bulk velocity. */
    int dimension() const override;

    /** This function returns true if the bulk velocity of the source is nonzero. */
    bool hasVelocity() const override;

    /** This function returns the wavelength range for this source. Outside this range, all
        luminosities are zero. This source's wavelength range is determined as the intersection of the
        simulation's source wavelength range (obtained from the simulation configuration) and the
        intrinsic wavelength range of the %SED associated with the source.

        This function implements the SourceWavelengthRangeInterface interface. */
    Range wavelengthRange() const override;

    /** This function returns the luminosity \f$L\f$ (i.e. radiative power) of the source
        integrated over the wavelength range of primary sources (configured for the source system
        as a whole) and across its complete spatial domain. */
     double luminosity() const override;

     /** This function returns the specific luminosity \f$L_\lambda\f$ (i.e. radiative power per
         unit of wavelength) of the source at the specified wavelength, or zero if the wavelength is
         outside the wavelength range of primary sources (configured for the source system as a
         whole) or if the source simply does not emit at the wavelength. */
     double specificLuminosity(double wavelength) const override;

     /** This function implements the BulkVelocityInterface interface. It returns the bulk velocity
         of this source, as configured by the user. */
     Vec bulkVelocity() const override;

     /** This function causes the photon packet \em pp to be launched from the source using the
         given history index and luminosity contribution. In this abstract class, the function
         handles the wavelength sampling and normalization, relying on the subclass to determine
         the position and propagation direction of the emission from the geometry of the source. */
    void launch(PhotonPacket* pp, size_t historyIndex, double L) const override;

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
    virtual void launchNormalized(PhotonPacket* pp, size_t historyIndex, double lambda, double Lw,
                                  BulkVelocityInterface* bvi) const = 0;

    //======================== Data Members ========================

private:
    // wavelength information initialized during setup
    bool _oligochromatic{false};    // true if the simulation is oligochromatic
    double _xi{0.};                 // the wavelength bias fraction
    WavelengthDistribution* _biasDistribution{nullptr}; // the wavelength bias distribution

    // pointer to an object offering the redshift interface; either "this" or null pointer if the bulk velocity is zero
    BulkVelocityInterface* _bvi{nullptr};
};

//////////////////////////////////////////////////////////////////////

#endif
