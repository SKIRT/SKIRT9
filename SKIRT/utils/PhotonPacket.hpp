/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PHOTONPACKET_HPP
#define PHOTONPACKET_HPP

#include "SpatialGridPath.hpp"
#include "StokesVector.hpp"
class AngularDistributionInterface;
class PolarizationProfileInterface;
class VelocityInterface;

////////////////////////////////////////////////////////////////////

/** The PhotonPacket class is used to describe monochromatic photon packets, the basic entities
    that are transported during a radiative transfer simulation. The fundamental properties of a
    photon packet include its \em wavelength and its \em weight. The wavelength property specifies
    the wavelength (or equivalenty, the frequency) of all photons in the packet. The weight
    property specifies the number of photons carried by the packet, or more precisely the number of
    photons per unit of time. Indeed, because we consider time-independent models, a photon packet
    actually represents a photon stream, consisting of a number of photons per unit of time. The
    weight can be fractional because it is derived from some arbitrary luminosity, and because it
    can be adjusted by arbitrary biasing factors during the photon packet's lifetime.

    Note that adjusting a photon packet's wavelength (perhaps because of a Doppler shift)
    indirectly affects the luminosity represented by the packet, because the latter is directly
    proportional to the frequency and thus inversely proportional to the wavelength.

    Apart from its wavelength and weight, a photon packet carries information about its polarization
    state, about its origin (e.g. emission by a primary or secondary source), about the
    interactions it experienced since its emission (e.g. the number of scattering events), and
    about its current path (e.g. starting position, propagation direction, spatial grid cells being
    crossed). A photon packet also carries a unique identifier (within the current emission
    segment) for the \em history it is part of, i.e. the collection of scattered and peeled-off
    packets derived from the same originally launched packet (emitted by a primary or secondary
    source). This history index can be used by instruments to group the fluxes of all photon
    packets belonging to the same history.

    Implementation notes
    --------------------

    For a given luminosity \f$L\f$ (representing the amount of energy carried by the photon packet
    per unit of time) and wavelength \f$\lambda=c/\nu\f$, the weight of the photon packet (i.e. the
    number of photons it carries per unit time) can be written as \f[ W = \frac{L}{h\nu} =
    \frac{L\lambda}{hc}. \f] To avoid multiplying and dividing by the constant factor \f$hc\f$,
    the weight tracked by a photon packet is instead simply defined as \f$W'=L\lambda\f$.

    The current path and the polarization state are handled by publicly inheriting repectively the
    SpatialGridPath class and the StokesVector class. For performance reasons, a PhotonPacket object
    is usually constructed once at the start of a loop and then reused in the loop body for many
    consecutive launches; this allows the vectors with path information to remain allocated. Also,
    several trivial functions are implemented inline in the header.

    The PhotonPacket class offers auxiliary facilities for caching "technical" information related
    to a photon packet's state during its lifecycle. This optional information is managed by the
    clients invoking the corresponding store and retrieve functions. Examples include the observed
    optical depth, used to avoid recalculating the optical depth for consecutive instruments with
    the same viewing direction, and extra scattering information, used by some material mixes to
    ensure that peel-off and random-walk operations for a given scattering interaction are treated
    consistently. */
class PhotonPacket : public SpatialGridPath, public StokesVector
{
    // ------- Construction, launch and lifecycle events -------

public:
    /** The constructor initializes an empty photon packet object. After construction, the photon
        packet is ready to be launched through one of the launch() functions. The other functions
        in this class should be invoked only after the photon packet has been launched. The same
        photon packet object can be re-launched multiple times. */
    PhotonPacket();

    /** This function initializes the photon packet for a new life cycle. The arguments specify the
        packet's history index, the wavelength of its photons, the source luminosity represented by
        the packet, its starting position, and its propagation direction. The function calculates
        the weight corresponding to the specified luminosity and wavelength, copies the remaining
        arguments to the corresponding data members and initializes the other data members as
        described below.

        If the VelocityInterface is specified (i.e. it is not the null pointer), then the
        packet's wavelength is Doppler shifted corresponding to its emission direction. If the
        PolarizationStateInterface is specified, the packet's polarization state is set according
        to its emission direction; otherwise the polarization state is set to unpolarized. The
        AngularDistributionInterface is not used by this function, but pointers to all three
        interfaces are stored in data members for use by the launchEmissionPeelOff() function.

        The emission origin is set to an invalid value, which should be overridden by calling the
        setPrimaryOrigin() or the setSecondaryOrigin() function. The number of scattering events is
        set to zero. The current path is invalidated, and all information about the previous life
        cycle is lost. */
    void launch(size_t historyIndex, double lambda, double L, Position bfr, Direction bfk,
                VelocityInterface* bvi = nullptr, AngularDistributionInterface* adi = nullptr,
                PolarizationProfileInterface* ppi = nullptr);

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
        derives and the direction towards the instrument. The function copies the relevant values
        from the base photon packet to the peel off photon packet, updates the peel off direction,
        and honors some extra source properties tracked by the base photon packet as follows.

        If the base packet has a VelocityInterface (i.e. if it is not the null pointer), then
        the peel-off packet's wavelength is Doppler shifted corresponding to its propagation
        direction. If the base packet has an AngularDistributionInterface, a bias for the
        probability of the peel-off propagation direction is applied to the peel-off packet's
        weight. If the base packet has a PolarizationStateInterface, the peel-off packet's
        polarization state is set according to its propagation direction; otherwise its
        polarization state is set to unpolarized.

        The current path of the peel off photon packet is invalidated, and all information about
        its previous life cycle is lost. The base photon packet remains unchanged. */
    void launchEmissionPeelOff(const PhotonPacket* pp, Direction bfk);

    /** This function initializes a peel off photon packet being sent to an instrument for a
        scattering event. The arguments specify the base photon packet from which the peel off
        derives, the direction towards the instrument, the new wavelength measured in a local frame
        with bulk velocity \f${\bf{v}}\f$, and the weight bias (as a multiplication factor). The
        function copies the relevant values from the base photon packet to the peel off photon
        packet, updates the peel off direction, the wavelength and weight, and increments the
        scattering counter.

        The peel off photon packet is initialized to an unpolarized state; the polarization state
        should be properly updated after the launch through the StokesVector class functions. The
        current path of the peel off photon packet is invalidated, and all information about its
        previous life cycle is lost. The base photon packet remains unchanged. */
    void launchScatteringPeelOff(const PhotonPacket* pp, Direction bfk, Vec bfv, double lambda, double w);

    /** This function causes the propagation of the photon packet over a physical distance \f$s\f$.
        It updates the position from \f${\bf{r}}\f$ to \f${\bf{r}}+s\,{\bf{k}}\f$, where
        \f${\bf{k}}\f$ is the propagation direction of the photon packet, invalidating the current
        path. */
    void propagate(double s);

    /** This function scatters the photon packet into the new direction \f${\bf{k}}\f$ with the new
        wavelength \f$\lambda\f$ measured in a local frame with bulk velocity \f${\bf{v}}\f$. It
        increments the counter that keeps track of scattering events, updates the direction,
        invalidating the current path, and updates the wavelength. The polarization remains
        unchanged; it should be properly updated through the StokesVector class functions. */
    void scatter(Direction bfk, Vec bfv, double lambda);

    /** This function applies the given weight bias given as a multiplication factor. */
    void applyBias(double w);

    // ------- Getting trivial properties -------

public:
    /** This function returns the wavelength \f$\lambda_0\f$ of the photon packet when it was
        launched, relative to the rest frame of the original source. */
    double sourceRestFrameWavelength() const { return _lambda0; }

    /** This function returns the current wavelength \f$\lambda\f$ of the photon packet relative to
        the model coordinate system. */
    double wavelength() const { return _lambda; }

    /** This function returns the luminosity \f$L\f$ represented by the photon packet, calculated
        from its current wavelength and weight. */
    double luminosity() const { return _W / _lambda; }

    /** This function returns true if the photon packet originated from a primary source, false
        otherwise. */
    bool hasPrimaryOrigin() const { return _compIndex > 0; }

    /** This function returns true if the photon packet originated from a medium, false otherwise.
        */
    bool hasSecondaryOrigin() const { return _compIndex < 0; }

    /** This function returns the index of the source or medium component that originally emitted
        the photon packet. */
    int compIndex() const { return abs(_compIndex) - 1; }

    /** This function returns the index of the source or medium component that originally emitted
        the photon packet. */
    size_t historyIndex() const { return _historyIndex; }

    /** This function returns the number of scattering events the photon packet has experienced. */
    int numScatt() const { return _nscatt; }

    // ------- Calculating Doppler shifts -------

public:
    /** This function returns the Doppler-shifted wavelength that should be assigned to a photon
        packet (i.e. the wavelength relative to the model coordinate frame) when the packet is
        emitted from a moving source (with non-relativistic velocity). The arguments specify the
        emitted wavelength \f$\lambda_\text{src}\f$ in the source rest frame, the direction
        \f${\bf{k}}_\text{ph}\f$ of the emitted photon packet, and the velocity of the source
        \f${\bf{v}}_\text{src}\f$ relative to the model coordinate frame. The photon packet
        wavelength \f$\lambda_\text{ph}\f$ can then be written as \f[ \lambda_\text{ph} =
        \lambda_\text{src} \left(1 - \frac{{\bf{k}}_\text{ph} \cdot {\bf{v}}_\text{src}}{c} \right)
        \f] where \f$c\f$ is the speed of light in vacuum. */
    static double shiftedEmissionWavelength(double sourceWavelength, Direction photonDirection, Vec sourceVelocity);

    /** This function returns the Doppler-shifted wavelength perceived by a moving receiver (with
        non-relativistic velocity) for an incoming photon packet. The arguments specify the photon
        packet wavelength \f$\lambda_\text{ph}\f$ in the model coordinate frame, the direction
        \f${\bf{k}}_\text{ph}\f$ of the incoming photon packet, the velocity of the receiver
        \f${\bf{v}}_\text{rec}\f$ relative to the model coordinate frame, and the expansion
        velocity \f$\Delta v_\mathrm{h}\f$ of the universe relative to the starting position of the
        photon packet's current path (or zero if cosmological expansion is not taken into account).

        The wavelength \f$\lambda_\text{rec}\f$ perceived by the receiver can then be written as
        \f[ \left. \lambda_\text{rec} = \lambda_\text{ph} \middle/ \left(1 -
        \frac{{\bf{k}}_\text{ph} \cdot {\bf{v}}_\text{rec} + \Delta v_\mathrm{h}}{c} \right)
        \right. \f] where \f$c\f$ is the speed of light in vacuum. */
    static double shiftedReceptionWavelength(double photonWavelength, Direction photonDirection, Vec receiverVelocity,
                                             double expansionVelocity = 0.);

    /** This function returns the Doppler-shifted wavelength perceived for this photon packet by a
        moving receiver (with non-relativistic velocity). The arguments specify the velocity of the
        receiver \f${\bf{v}}_\text{rec}\f$ relative to the model coordinate frame and the expansion
        velocity \f$\Delta v_\mathrm{h}\f$ of the universe relative to the starting position of the
        photon packet's current path (or zero if cosmological expansion is not taken into account).
        See the shiftedReceptionWavelength() function for more details. */
    double perceivedWavelength(Vec receiverVelocity, double expansionVelocity) const;

    /** This function returns the luminosity \f$L\f$ represented by the photon packet, calculated
        from its current weight and the specified perceived wavelength. */
    double perceivedLuminosity(double lambda) const { return _W / lambda; }

    // ------- Caching observed optical depth -------

public:
    /** This function stores the most recently "observed" optical depth, calculated externally, in
        a data member. This capability is offered so that consecutive instruments with the same
        observer type, position and viewing direction can avoid recalculating the optical depth. */
    void setObservedOpticalDepth(double tau)
    {
        _observedOpticalDepth = tau;
        _hasObservedOpticalDepth = true;
    }

    /** This function returns true if an "observed" optical depth value has been stored and the
        packet has not since been relaunched. Otherwise the function returns false. */
    bool hasObservedOpticalDepth() const { return _hasObservedOpticalDepth; }

    /** If hasObservedOpticalDepth() returns true, this function returns the most recently stored
        "observed" optical depth. Otherwise, it returns some meaningless value. This capability is
        offered so that consecutive instruments with the same observer type, position and viewing
        direction can avoid recalculating the optical depth. */
    double observedOpticalDepth() const { return _observedOpticalDepth; }

    // ------- Caching scattering info -------

public:
    /** This function sets the index of the medium component participating in the next peel-off or
        random-walk scattering operation. It should be called by the photon cycle machinery just
        before the material mix corresponding to the specified medium component is requested to
        determine its contribution to a peel-off or to perform a random-walk scattering operation.
        This mechanism allows the photon packet to maintain a seperate scattering information
        record for each medium component. */
    void setScatteringComponent(int h) { _h = h; }

    /** This class defines the scattering information record maintained by a photon packet for each
        medium component that needs it. In case some aspects of the scattering event are governed
        by a probability distribution, this information record can be used to ensure that the
        outgoing random walk photon packet and all corresponding peel-off photon packets use the
        same randomly sampled values.

        The following data members have public scope so that they can be directly accessed by
        client code:

        - \em valid: boolean flag, initially false; should be set to true once this record holds
        valid values.

        - \em dipole: boolean flag, initially false, indicating the selected scattering phase
        function: true means dipole, false means isotropic.

        - \em species: integer identifier of the scattering particle species, initially zero.

        - \em velocity: velocity vector of the scattering particle in the local frame, initially
        the null vector.

        - \em lambda: outgoing wavelength of the scattered photon packet, initially zero.

        Apart from initialization by the getScatteringInfo() function as indicated above, the
        client code is fully responsible for writing and reading these fields. */
    class ScatteringInfo
    {
    private:
        friend class PhotonPacket;
        int _h;  // the index of the medium component to which this information record belongs

    public:
        ScatteringInfo(int h) : _h{h} {}  // cannot be private because we use emplace_back to construct

    public:
        bool valid{false};   // true if this record holds valid values
        bool dipole{false};  // true if scattering as a dipole, false if scattering isotropically
        int species{0};      // identifier of the scattering particle species
        Vec velocity;        // the velocity vector of the scattering particle in the local frame
        double lambda{0.};   // the outgoing wavelength of the scattered photon packet
    };

    /** This function returns a pointer to the (writable) scattering information record maintained
        for the current medium component as previously set through the setScatteringComponent()
        function. When this function is called for the first time for a new scattering event, the
        information record will be freshly initialized as described for the ScatteringInfo class.
        Specifically, the \em valid flag will be false. Subsequent invocations for the same
        scattering event will return the same information record with field values as updated by
        the client code.

        Material mixes that need to store scattering information should call this function for each
        peel-off and random-walk scattering operation. If \em valid is false, the material mix
        should update the relevant fields and set \em valid to true. If \em valid is already true,
        the previously stored field values can be retrieved. */
    ScatteringInfo* getScatteringInfo();

    // ------- Data members -------

private:
    // current physical properties (in addition to inherited data members)
    double _lambda{0};  // current wavelength relative to the model coordinate system
    double _W{0};       // current weight, defined as L*lambda to avoid division and multiplication by hc

    // physical information on radiation source; the interfaces are not used in peel-off photon packets
    double _lambda0{0};  // original wavelength in the rest-frame of the source
    VelocityInterface* _bvi{nullptr};
    AngularDistributionInterface* _adi{nullptr};
    PolarizationProfileInterface* _ppi{nullptr};

    // administrative information on origin
    int _compIndex{0};        // sign * (index of the originating source or medium component + 1)
                              //  0: uninitialized   >0: primary   <0: secondary
    size_t _historyIndex{0};  // index of the photon packet's history in the current emission segment

    // information on life cycle
    int _nscatt{0};  // number of experienced scattering events

    // observed optical depth
    double _observedOpticalDepth{0.};      // optical depth calculated for peel-off to an instrument
    bool _hasObservedOpticalDepth{false};  // true if the above field holds a valid value for this packet

    // scattering information
    int _h{0};  // the index of the medium component currently participating in a scattering operation
    vector<ScatteringInfo> _scatteringInfo;  // list of scattering information records
};

////////////////////////////////////////////////////////////////////

#endif
