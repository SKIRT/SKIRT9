/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SOURCE_HPP
#define SOURCE_HPP

#include "SimulationItem.hpp"
#include "SourceWavelengthRangeInterface.hpp"
#include "WavelengthDistribution.hpp"
class PhotonPacket;
class Random;

//////////////////////////////////////////////////////////////////////

/** Source is an abstract class representing a primary radiation source in the simulation.
    Each source (an instance of a Source subclass) must define the following information at each
    point in the spatial domain of the source (which is independent of the spatial grid in the
    simulation):
      - The spectral energy distribution (SED) of the emission averaged over the unit sphere.
      - Some normalization for the luminosity (e.g., bolometric, at a given wavelength, ...).
      - If the emission is anisotropic, the rest-frame angular distribution of the emission.
      - If the emission is polarized, the polarization state of the emission in each direction.
      - The velocity of the source relative to the model coordinate frame.

    Furthermore, each source has a function for launching a photon packet that proceeds roughly
    as follows:
      - Sample a location from the spatial density distribution.
      - Sample a wavelength from the SED at that location (also see below).
      - Determine the corresponding (possibly biased) luminosity weight
      - Determine the rest-frame angular distribution of the emission at that location and at
        that wavelength, i.e. an object offering the AngularDistribution interface (functions
        to return a random direction, and to return the probability for a given direction)
      - Determine the polarization state of the emission at that location and at that wavelength,
        i.e. an object offering the PolarizationState interface (function to return a Stokes vector
        for a given direction)
      - Determine the velocity of the source at that location
      - Pass the items listed above to the photon packet launch procedure

    Wavelengths for new photon packets can be sampled from the intrinsic spectral distribution of
    the source \f$s(\lambda)\f$ and/or from a \em bias wavelength distribution \f$b(\lambda)\f$. Both
    the bias fraction \f$\xi\f$ and the bias distribution \f$b(\lambda)\f$ can be configured by the
    user. Given these distributions and bias factor, the composite distribution \f$q(\lambda)\f$ is

    \f[ q(\lambda) = (1-\xi) s(\lambda) + \xi b(\lambda) \f]

    and the corresponding biasing weight factor becomes

    \f[ w(\lambda) = \frac{s(\lambda)}{q(\lambda)} = \frac{s(\lambda)}{(1-\xi) s(\lambda) + \xi b(\lambda)} \f]

    By default, half of the photon packet wavelengths are sampled from each of the distributions,
    and the default bias distribution spreads wavelengths logarithmically over the wavelength range
    of the source (more precisely, the logarithm of the wavelength is distributed uniformly).
*/
class Source : public SimulationItem, public SourceWavelengthRangeInterface
{
    ITEM_ABSTRACT(Source, SimulationItem, "a primary radiation source")

        ATTRIBUTE_SUB_PROPERTIES_HERE(Source)

        PROPERTY_DOUBLE(sourceWeight, "the weight of this source for the number of photon packets launched")
        ATTRIBUTE_MIN_VALUE(sourceWeight, "]0")
        ATTRIBUTE_MAX_VALUE(sourceWeight, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(sourceWeight, "1")
        ATTRIBUTE_DISPLAYED_IF(sourceWeight, "Level3")

        PROPERTY_DOUBLE(wavelengthBias, "the fraction of photon packet wavelengths sampled from a bias distribution")
        ATTRIBUTE_MIN_VALUE(wavelengthBias, "[0")
        ATTRIBUTE_MAX_VALUE(wavelengthBias, "1]")
        ATTRIBUTE_DEFAULT_VALUE(wavelengthBias, "0.5")
        ATTRIBUTE_RELEVANT_IF(wavelengthBias, "Panchromatic")
        ATTRIBUTE_DISPLAYED_IF(wavelengthBias, "Level3")

        PROPERTY_ITEM(wavelengthBiasDistribution, WavelengthDistribution,
                      "the bias distribution for sampling photon packet wavelengths")
        ATTRIBUTE_DEFAULT_VALUE(wavelengthBiasDistribution, "LogWavelengthDistribution")
        ATTRIBUTE_RELEVANT_IF(wavelengthBiasDistribution, "Panchromatic&wavelengthBias")
        ATTRIBUTE_DISPLAYED_IF(wavelengthBiasDistribution, "Level3")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function caches the simulation's random generator for use by subclasses. */
    void setupSelfBefore() override;

    /** This function logs a warning message if the given range is smaller than the configured
        source wavelength range. The second argument specifies the type of the simulation item to
        be included in the message. This function can be called from subclasses. */
    void informAvailableWavelengthRange(Range available, string itemType);

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the source, which depends on the (lack of) symmetry
        of its geometry. A value of 1 means spherical symmetry, 2 means axial symmetry and 3 means
        none of these symmetries. */
    virtual int dimension() const = 0;

    /** This function returns true if this source may have a nonzero velocity for some positions.
        It may be called before setup of the receiving source has completed. */
    virtual bool hasVelocity() const = 0;

    /** This function returns the luminosity \f$L\f$ (i.e. radiative power) of the source
        integrated over the wavelength range of primary sources (configured for the source system
        as a whole) and across its complete spatial domain. */
    virtual double luminosity() const = 0;

    /** This function returns the specific luminosity \f$L_\lambda\f$ (i.e. radiative power per
        unit of wavelength) of the source at the specified wavelength, or zero if the wavelength is
        outside the wavelength range of primary sources (configured for the source system as a
        whole) or if the source simply does not emit at the wavelength. */
    virtual double specificLuminosity(double wavelength) const = 0;

    /** This function may perform some preparations for launching photon packets. It is called in
        serial mode before each segment of photon packet launches, providing the history indices
        mapped by the source system to this particular source. This allows the source to further
        map these indices to subsources if that would be meaningful. See the description of the
        SourceSystem class for more information. The default implementation of this function does
        nothing. */
    virtual void prepareForLaunch(double sourceBias, size_t firstIndex, size_t numIndices);

    /** This function causes the photon packet \em pp to be launched from the source using the
        given history index and luminosity contribution. The photon packet's contents is fully
        (re-)initialized so that it is ready to start its lifecycle. */
    virtual void launch(PhotonPacket* pp, size_t historyIndex, double L) const = 0;

    //======================== Other Functions =======================

protected:
    /** This function returns the simulation's random generator as a service to subclasses. */
    Random* random() const { return _random; }

    //======================== Data Members ========================

private:
    // data members initialized during setup by setupSelfBefore()
    Random* _random{nullptr};
};

//////////////////////////////////////////////////////////////////////

#endif
