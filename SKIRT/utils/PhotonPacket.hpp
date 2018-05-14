/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PHOTONPACKET_HPP
#define PHOTONPACKET_HPP

#include "DustGridPath.hpp"
#include "StokesVector.hpp"

////////////////////////////////////////////////////////////////////

/** The PhotonPacket class is used to describe photon packets, the basic luminosity packets that
    are transported during a radiative transfer simulation. Photon packets are monochromatic, i.e.
    they contain photons at a single wavelength. Apart from its luminosity, wavelength and
    polarization state, a photon packet carries information about its origin (e.g. emission by a
    primary or secondary source), about the interactions it experienced since its emission (e.g.
    the number of scattering events), and about its current path (e.g. starting position,
    propagation direction, list of dust cells being crossed). A photon packet also carries a unique
    identifier (within the current emission segment) for the \em history it is part of, i.e. the
    collection of scattered and peeled-off packets derived from the same originally launched packet
    (emitted by a primary or secondary source). This history index is used by some instruments to
    group the fluxes of all photon packets belonging to the same history.

    Because a photon packet's wavelength and/or frequency are retrieved much more frequently than
    updated, the PhotonPacket class has a data member for both both wavelength and frequency. The
    setters for wavelength and frequency immediately calculate and store both values, so that the
    getters can trivially return the stored value.

    Furthermore, a photon packet carries a specific luminosity rather than a luminosity. In line
    with common practice, we use \f$L_\nu(\lambda)\f$ expressed in units of W/Hz. An important
    benefit of this convention is that, when updating the frequency for a Doppler shift, the
    specific luminosity remains constant, because the frequency interval scales with the frequency.

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
        packet's history index, the specific luminosity, the wavelength, the starting position and
        the propagation direction. The function copies these values to the corresponding data
        members and initializes the other data members to default values. The number of scattering
        events is set to zero. The emission origin is set to an invalid value, which should be
        overridden by calling the setPrimaryOrigin() or the setSecondaryOrigin() function. The
        polarization state is set to unpolarized, which can be overridden through the
        setPolarized() function. The current path is invalidated by these changes, and all
        information about the previous life cycle is lost. */
    void launch(size_t historyIndex, double Lnu, double lambda, Position bfr, Direction bfk);

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
        derives, the direction towards the instrument, and the luminosity bias in case of
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
        derives, the direction towards the instrument, and the luminosity bias (as a multiplication
        factor). The function copies the relevant values from the base photon packet to the peel
        off photon packet, updates the peel off direction and luminosity, and increments the
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

    /** This function applies the luminosity bias given as a multiplication factor. */
    void applyLuminosityBias(double w);

    /** This function adjusts the packet's frequency and wavelength according to the specified
        (non-relativistic) redshift \f$z=v/c\f$: \f[\begin{aligned} \lambda' &= \lambda (1+z) \\
        \nu' &= \nu / (1+z) \end{aligned}\f] In other words, if \f$z>0\f$, then the packet is
        shifted to longer wavelengths (and lower frequencies). */
    void applyRedshift(double z);


    // ------- Getting trivial properties -------

    /** This function returns the specific luminosity \f$L_\nu(\lambda)\f$ of the photon packet. */
    double luminosity() const { return _Lnu; }

    /** This function returns the wavelength \f$\lambda\f$ of the photon packet. */
    double wavelength() const { return _lambda; }

    /** This function returns the wavelength \f$\lambda\f$ of the photon packet. */
    double lambda() const { return _lambda; }

    /** This function returns the frequency \f$\nu\f$ of the photon packet. */
    double frequency() const { return _nu; }

    /** This function returns the frequency \f$\nu\f$ of the photon packet. */
    double nu() const { return _nu; }

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
    double _Lnu{0};          // specific luminosity (W/Hz)
    double _lambda{0};       // wavelength (m)
    double _nu{0};           // frequency (Hz)

    // information on history
    int _nscatt{0};          // number of experienced scattering events

    // information on origin
    int _compIndex{0};       // sign * (index of the originating source or medium component + 1)
                             //  0: uninitialized   >0: primary   <0: secondary
    size_t _historyIndex{0}; // index of the photon packet's history in the current emission segment
};

////////////////////////////////////////////////////////////////////

#endif
