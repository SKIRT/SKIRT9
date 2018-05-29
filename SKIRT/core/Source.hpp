/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SOURCE_HPP
#define SOURCE_HPP

#include "SimulationItem.hpp"
class PhotonPacket;

//////////////////////////////////////////////////////////////////////

/** Source is an abstract class representing a primary radiation source in the simulation.
    Each source (an instance of a Source subclass) must define the following information at each
    point in the spatial domain of the source (which is independent of the spatial grid in the
    simulation):
      - The spectral energy distribution (SED) of the emission averaged over the unit sphere.
      - Some normalization for the luminosity (e.g., bolometric, at a given wavelength, ...).
      - If the emission is anisotropic, the rest-frame angular distribution of the emission.
      - If the emission is polarized, the polarization state of the emission in each direction.
      - The bulk velocity of the source relative to the model coordinate frame.

    Furthermore, each source has a function for launching a photon packet that proceeds roughly
    as follows:
      - Sample a location from the spatial density distribution.
      - Sample a wavelength from the SED at that location.
      - Determine the corresponding (possibly biased) luminosity weight
      - Determine the rest-frame angular distribution of the emission at that location and at
        that wavelength, i.e. an object offering the AngularDistribution interface (functions
        to return a random direction, and to return the probability for a given direction)
      - Determine the polarization state of the emission at that location and at that wavelength,
        i.e. an object offering the PolarizationState interface (function to return a Stokes vector
        given a propagation direction)
      - Determine the bulk velocity of the source at that location
      - Pass the items listed above to the photon packet launch procedure
*/
class Source : public SimulationItem
{
    ITEM_ABSTRACT(Source, SimulationItem, "a primary radiation source")

    PROPERTY_DOUBLE(emissionWeight, "the weight of this source for the number of photon packets launched")
        ATTRIBUTE_MIN_VALUE(emissionWeight, "]0")
        ATTRIBUTE_MAX_VALUE(emissionWeight, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(emissionWeight, "1")
        ATTRIBUTE_SILENT(emissionWeight)

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the source, which depends on the (lack of) symmetry
        of its geometry. A value of 1 means spherical symmetry, 2 means axial symmetry and 3 means
        none of these symmetries. */
    virtual int dimension() const = 0;

    /** This function returns the bolometric luminosity \f$L\f$ of the source across its spatial
        and spectral domain. */
    virtual double luminosity() const = 0;

    /** This function may perform some preparations for launching photon packets. It is called in
        serial mode before each segment of photon packet launches, providing the history indices
        mapped by the source system to this particular source. This allows the source to further
        map these indices to subsources if that would be meaningful. See the description of the
        SourceSystem class for more information. The default implementation of this function does
        nothing. */
    virtual void prepareForLaunch(size_t firstIndex, size_t numIndices);

    /** This function causes the photon packet \em pp to be launched from the source using the
        given history index and luminosity contribution. The photon packet's contents is fully
        (re-)initialized so that it is ready to start its lifecycle. */
    virtual void launch(PhotonPacket* pp, size_t historyIndex, double L) const = 0;
};

//////////////////////////////////////////////////////////////////////

#endif
