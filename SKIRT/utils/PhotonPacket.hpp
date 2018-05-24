/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PHOTONPACKET_HPP
#define PHOTONPACKET_HPP

#include "DustGridPath.hpp"
#include "StokesVector.hpp"

////////////////////////////////////////////////////////////////////

/** The PhotonPacket class is used to describe monochromatic photon packets, the basic entities
    that are transported during a radiative transfer simulation. The fundamental properties of a
    photon packet include its \em frequency and its \em weight. The frequency property specifies
    the frequency (or equivalenty, the wavelength) of all photons in the packet. The weight
    property specifies the number of photons carried by the packet, or more precisely the number of
    photons per unit of time. Indeed, because we consider time-independent models, a photon packet
    actually represents a photon stream, consisting of a number of photons per unit of time. The
    weight can be fractional because it is derived from some arbitrary luminosity, and because it
    can be adjusted by arbitrary biasing factors during the photon packet's lifetime.

    Apart from its frequency and weight, a photon packet carries information about its polarization
    state, about its origin (e.g. emission by a primary or secondary source), about the
    interactions it experienced since its emission (e.g. the number of scattering events), and
    about its current path (e.g. starting position, propagation direction, list of dust cells being
    crossed). A photon packet also carries a unique identifier (within the current emission
    segment) for the \em history it is part of, i.e. the collection of scattered and peeled-off
    packets derived from the same originally launched packet (emitted by a primary or secondary
    source). This history index can be used by instruments to group the fluxes of all photon
    packets belonging to the same history.

    Because a photon packet's wavelength and/or frequency are retrieved much more frequently than
    updated, the PhotonPacket class has a data member for both both wavelength and frequency. The
    setters for wavelength and frequency immediately calculate and store both values, so that the
    getters can trivially return the stored value.

    The current path and the polarization state are handled by publicly inheriting repectively the
    DustGridPath class and the StokesVector class. For performance reasons, a PhotonPacket object
    is usually constructed once at the start of a loop and then reused in the loop body for many
    consecutive launches; this allows the vectors with path information to remain allocated. Also,
    some trivial functions are implemented inline in the header. */
class PhotonPacket : public DustGridPath, public StokesVector
{
public:

    // ------- Construction, launch and lifecycle events -------

    /** The constructor initializes an empty photon packet object. After construction, the photon
        packet is ready to be launched through one of the launch() functions. The other functions
        in this class should be invoked only after the photon packet has been launched. The same
        photon packet object can be re-launched multiple times. */
    PhotonPacket();

    /** This function initializes the photon packet for a new life cycle. The arguments specify the
        packet's history index, the wavelength of its photons, the source luminosity represented by
        the packet, its starting position, and its propagation direction. The function calculates
        the weight corresponding to the specified luminosity and wavelength, copies the remaining
        arguments to the corresponding data members and initializes the other data members to
        default values. The number of scattering events is set to zero. The emission origin is set
        to an invalid value, which should be overridden by calling the setPrimaryOrigin() or the
        setSecondaryOrigin() function. The polarization state is set to unpolarized, which can be
        overridden through the setPolarized() function. The current path is invalidated by these
        changes, and all information about the previous life cycle is lost. */
    void launch(size_t historyIndex, double lambda, double L, Position bfr, Direction bfk);

    /** This function establishes that the photon packet has been emitted by a primary source and
        registers the index of the emitting source component. This information is used by some
        instruments to record fluxes seperately based on their origin. This function should be
        called just after launch. */
    void setPrimaryOrigin(int sourceCompIndex);

    /** This function establishes that the photon packet has been emitted by a medium and registers
        the index of the emitting medium component. This information is used by some instruments to
        record fluxes seperately based on their origin. This function should be called just after
        launch. */
    void setSecondaryOrigin(int mediumCompIndex);

    /** This function initializes a peel off photon packet being sent to an instrument for an
        emission event. The arguments specify the base photon packet from which the peel off
        derives, the direction towards the instrument, and the weight bias in case of
        anisotropic emission (as a multiplication factor). The function copies the relevant values
        from the base photon packet to the peel off photon packet, updates the peel off direction,
        and applies the bias. The peel off photon packet is initialized to an unpolarized state;
        the polarization state should be properly updated after the launch through the StokesVector
        class functions. The current path of the peel off photon packet is invalidated, and all
        information about its previous life cycle is lost. The base photon packet remains
        unchanged. */
    void launchEmissionPeelOff(const PhotonPacket* pp, Direction bfk, double w);

    /** This function initializes a peel off photon packet being sent to an instrument for a
        scattering event. The arguments specify the base photon packet from which the peel off
        derives, the direction towards the instrument, and the weight bias (as a multiplication
        factor). The function copies the relevant values from the base photon packet to the peel
        off photon packet, updates the peel off direction and weight, and increments the
        scattering counter. The peel off photon packet is initialized to an unpolarized state; the
        polarization state should be properly updated after the launch through the StokesVector
        class functions. The current path of the peel off photon packet is invalidated, and all
        information about its previous life cycle is lost. The base photon packet remains
        unchanged. */
    void launchScatteringPeelOff(const PhotonPacket* pp, Direction bfk, double w);

    /** This function causes the propagation of the photon packet over a physical distance \f$s\f$.
        It updates the position from \f${\bf{r}}\f$ to \f${\bf{r}}+s\,{\bf{k}}\f$, where
        \f${\bf{k}}\f$ is the propagation direction of the photon packet, invalidating the current
        path. */
    void propagate(double s);

    /** This function scatters the photon packet into the new direction \f${\bf{k}}\f$. It
        increments the counter that keeps track of scattering events and updates the direction,
        invalidating the current path. The polarization remains unchanged; it should be properly
        updated through the StokesVector class functions. */
    void scatter(Direction bfk);

    /** This function applies the weight bias given as a multiplication factor. */
    void applyBias(double w);

    /** This function adjusts the packet's frequency and wavelength according to the specified
        (non-relativistic) redshift \f$z=v/c\f$: \f[\begin{aligned} \lambda' &= \lambda (1+z) \\
        \nu' &= \nu / (1+z) \end{aligned}\f] In other words, if \f$z>0\f$, then the packet is
        shifted to longer wavelengths (and lower frequencies).

        Note that calling this function with nonzero redshift indirectly affects the luminosity
        represented by the packet, because the latter is directly proportional to the frequency. */
    void applyRedshift(double z);

    // ------- Getting trivial properties -------

    /** This function returns the wavelength \f$\lambda\f$ of the photon packet. */
    double wavelength() const { return _lambda; }

    /** This function returns the wavelength \f$\lambda\f$ of the photon packet. */
    double lambda() const { return _lambda; }

    /** This function returns the frequency \f$\nu\f$ of the photon packet. */
    double frequency() const { return _nu; }

    /** This function returns the frequency \f$\nu\f$ of the photon packet. */
    double nu() const { return _nu; }

    /** This function returns the luminosity \f$L\f$ represented by the photon packet, calculated
        from its current frequency and weight. */
    double luminosity() const { return _W * _nu; }

    /** This function returns the number of scattering events the photon packet has experienced. */
    int numScatt() const { return _nscatt; }

    /** This function returns true if the photon packet originated from a primary source, false
        otherwise. */
    bool hasPrimaryOrigin() const { return _compIndex > 0; }

    /** This function returns true if the photon packet originated from a medium, false otherwise.
        */
    bool hasSecondaryOrigin() const { return _compIndex < 0; }

    /** This function returns the index of the source or medium component that originally emitted
        the photon packet. */
    int compIndex() const { return abs(_compIndex)-1; }

    /** This function returns the index of the source or medium component that originally emitted
        the photon packet. */
    size_t historyIndex() const { return _historyIndex; }

    // ------- Data members -------

private:
    // physical properties
    //  - the wavelength and frequency are relative to the model coordinate system
    //  - both values are always kept in sync by the corresponding setters
    double _lambda{0};       // wavelength (m)
    double _nu{0};           // frequency (Hz)
    double _W{0};            // weight (J), defined as L/nu to avoid division and multiplication by Planck constant

    // information on history
    int _nscatt{0};          // number of experienced scattering events

    // information on origin
    int _compIndex{0};       // sign * (index of the originating source or medium component + 1)
                             //  0: uninitialized   >0: primary   <0: secondary
    size_t _historyIndex{0}; // index of the photon packet's history in the current emission segment
};

////////////////////////////////////////////////////////////////////

#endif
