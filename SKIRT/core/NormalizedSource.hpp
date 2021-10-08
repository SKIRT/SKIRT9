/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NORMALIZEDSOURCE_HPP
#define NORMALIZEDSOURCE_HPP

#include "LuminosityNormalization.hpp"
#include "SED.hpp"
#include "Source.hpp"
class ContSED;

//////////////////////////////////////////////////////////////////////

/** NormalizedSource is an abstract class representing a primary radiation source characterized by
    a single SED object, i.e. the spectral distribution is identical in all spatial locations. The
    bolometric power of the source is characterized by a LuminosityNormalization object.

    Subclasses must handle the spatial distribution of the source and can optionally add bulk
    velocity, anisotropy and/or polarization. */
class NormalizedSource : public Source
{
    ITEM_ABSTRACT(NormalizedSource, Source, "a primary source with a single SED")

        ATTRIBUTE_SUB_PROPERTIES_HERE(NormalizedSource)

        PROPERTY_ITEM(sed, SED, "the spectral energy distribution for the source")
        ATTRIBUTE_DEFAULT_VALUE(sed, "BlackBodySED")

        PROPERTY_ITEM(normalization, LuminosityNormalization, "the type of luminosity normalization for the source")
        ATTRIBUTE_DEFAULT_VALUE(normalization, "IntegratedLuminosityNormalization")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function caches some wavelength information. */
    void setupSelfBefore() override;

    /** This function warns the user if this source's intrinsic wavelength range does not fully
        cover the configured wavelength range. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
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
        whole) or if the source simply does not emit at the wavelength. Discrete line emission
        (i.e. with zero line width) is ignored. */
    double specificLuminosity(double wavelength) const override;

    /** This function causes the photon packet \em pp to be launched from the source using the
         given history index and luminosity contribution. In this abstract class, the function
         handles the wavelength sampling and normalization, relying on the subclass to determine
         the position and propagation direction of the emission from the geometry of the source. */
    void launch(PhotonPacket* pp, size_t historyIndex, double L) const override;

    //============== Functions to be implemented in each subclass =============

    /** This function causes the photon packet \em pp to be launched from the source using the
        given history index, wavelength and weighted luminosity contribution. It must be
        implemented in a subclass to handle the spatial distribution of the source, optionally
        adding bulk velocity, anisotropy and/or polarization. */
    virtual void launchNormalized(PhotonPacket* pp, size_t historyIndex, double lambda, double Lw) const = 0;

    //======================== Data Members ========================

private:
    // wavelength information initialized during setup
    bool _oligochromatic{false};                         // true if the simulation is oligochromatic
    double _xi{0.};                                      // the wavelength bias fraction
    WavelengthDistribution* _biasDistribution{nullptr};  // the wavelength bias distribution
    ContSED* _contsed{nullptr};                          // cast of SED to ContSED or nullptr if line spectrum
};

//////////////////////////////////////////////////////////////////////

#endif
