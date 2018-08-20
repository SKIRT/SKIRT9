/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MONTECARLOSIMULATION_HPP
#define MONTECARLOSIMULATION_HPP

#include "Simulation.hpp"
#include "Configuration.hpp"
#include "InstrumentSystem.hpp"
#include "MediumSystem.hpp"
#include "ProbeSystem.hpp"
#include "SimulationMode.hpp"
#include "SourceSystem.hpp"
#include <atomic>

//////////////////////////////////////////////////////////////////////

/** The MonteCarloSimulation class is the top-level class describing a SKIRT simulation. Running a
    Monte Carlo simulation with SKIRT essentially comes down to constructing an instance of the
    MonteCarloSimulation class and invoking the setupAndRun() function on it.

    The MonteCarloSimulation class holds the source, media, instrument and probe systems,
    implements the core aspects of the photon packet life-cycle, and manages the iterative
    processes in the simulation (including phases, iterations and segments).

    TODO: more documentation ...

    The MonteCarloSimulation class also holds the non-discoverable \em config property, which is
    automatically set to an instance of the Configuration class. The setup() function of the config
    object is invoked at the very early stages of overall simulation setup, so that it can
    initialize its internal state to reflect the simulation configuration. As a result, it is safe
    for other simulation items to retrieve information from the config object during setup. */
class MonteCarloSimulation : public Simulation
{
    ITEM_CONCRETE(MonteCarloSimulation, Simulation, "a Monte Carlo simulation")

    PROPERTY_ITEM(mode, SimulationMode, "the overall simulation mode")
        ATTRIBUTE_DEFAULT_VALUE(mode, "ExtinctionOnlyMode")

    PROPERTY_ITEM(sourceSystem, SourceSystem, "the source system")
        ATTRIBUTE_DEFAULT_VALUE(sourceSystem, "SourceSystem")

    PROPERTY_ITEM(mediumSystem, MediumSystem, "the medium system")
        ATTRIBUTE_DEFAULT_VALUE(mediumSystem, "MediumSystem")
        ATTRIBUTE_OPTIONAL(mediumSystem)

    PROPERTY_ITEM(instrumentSystem, InstrumentSystem, "the instrument system")
        ATTRIBUTE_DEFAULT_VALUE(instrumentSystem, "InstrumentSystem")

    PROPERTY_ITEM(probeSystem, ProbeSystem, "the probe system")
        ATTRIBUTE_DEFAULT_VALUE(probeSystem, "ProbeSystem")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function performs setup for the complete simulation hierarchy. It calls the regular
        setup() function and notifies the probe system when setup has been completed. */
    void setupSimulation() override;

    /** This function performs initial setup for the MonteCarloSimulation object; it caches some
        frequently used pointers. */
    void setupSelfBefore() override;

    /** This function performs final setup for the MonteCarloSimulation object; it logs the
        dimension of the simulation. */
    void setupSelfAfter() override;

    //======== Getters for Non-Discoverable Properties =======

public:
    /** Returns the Configuration object for this simulation hierarchy. */
    Configuration* config() const;

    //======================== Other Functions =======================

protected:
    /** This function actually runs the simulation, assuming setup has been completed. */
    void runSimulation() override;

private:
    /** In a multi-processing environment, this function logs a message and waits for all processes
        to finish the work (i.e. it places a barrier). The string argument is included in the log
        message to indicate the scope of work that is being finished. If there is only a single
        process, the function does nothing. */
    void wait(string scope);

    /** This function initializes the progress counter used in logprogress() for the specified
        segment and logs the number of photon packets to be processed. */
    void initProgress(string segment, size_t numTotal);

    /** This function logs a progress message for the segment specified in the initprogress()
        function if the previous message was issued at least 3 seconds ago. The function
        must be called regularly while processing photon packets. The argument specifies the
        number of photon packets processed so far. */
    void logProgress(size_t numDone);

    /** This function launches the specified chunk of photon packets from primary sources. */
    void doPrimaryEmissionChunk(size_t firstIndex, size_t numIndices);

    /** This function simulates the peel-off of a photon packet after an emission event. This
        means that we create peel-off or shadow photon packets, one for every instrument in the
        instrument system, that we force to propagate in the direction of the observer(s) instead
        of in the propagation direction \f${\bf{k}}\f$ determined randomly by the emission process.
        Each peel-off photon packet has the same characteristics as the original peel-off photon
        packet, except that the propagation direction is altered to the direction
        \f${\bf{k}}_{\text{obs}}\f$ of the observer. For anistropic emission, a weight factor is
        applied to the luminosity to compensate for the fact that the probability that a photon
        packet would have been emitted towards the observer is not the same as the probability
        that it is emitted in any other direction. For each instrument in the instrument system,
        the function creates such a peel-off photon packet and feeds it to the instrument. The
        first argument specifies the photon packet that was just emitted; the second argument
        provides a placeholder peel off photon packet for use by the function. */
    void peelOffEmission(const PhotonPacket* pp, PhotonPacket* ppp);

    /** This function simulates the escape from the system and the absorption by dust of a fraction
        of the luminosity of a photon packet. It actually splits the luminosity \f$L_\ell\f$ of
        the photon packet in \f$N+2\f$ different parts, with \f$N\f$ the number of dust cells
        along its path through the dust system: a part \f$L_\ell^{\text{esc}}\f$ that escapes from
        the system, a part \f$L_\ell^{\text{sca}}\f$ that is scattered before it leaves the system,
        and \f$N\f$ parts \f$L_{\ell,n}^{\text{abs}}\f$, each of them corresponding to the fraction
        that is absorbed in a dust cell along the path. The part that scatters is the actual part
        of the photon packet that survives and continues in the photon packet's life cycle. So we
        alter the luminosity of the photon packet from \f$L_\ell\f$ to \f$L_\ell^{\text{sca}}\f$.
        If the absorption rates are not considered in the simulation, that is also the only part
        that matters.

        On the other hand, if the absorption rates need to be taken into account in the simulation
        (if dust emission is included in the simulation, or the mean radiation field is
        calculated), we have to simulate the absorption of the fractions
        \f$L_{\ell,n}^{\text{abs}}\f$ by the appropriate dust cells. It is obviously crucial to
        have expressions for \f$L_\ell^{\text{sca}}\f$ and \f$L_{\ell,n}^{\text{abs}}\f$ (and
        \f$L_\ell^{\text{esc}}\f$, but this does not really matter...). If we denote the total
        optical depth along the path of the photon packet as \f$\tau_{\ell,\text{path}}\f$ (this
        quantity is stored in the PhotonPacket object provided as an input parameter of this
        function), the fraction of the luminosity that escapes from the system without any
        interaction is \f[ L_\ell^{\text{esc}} = L_\ell\, {\text{e}}^{-\tau_{\ell,\text{path}}}.
        \f] Remains to subdivide the remainder of the initial luminosity, \f[ L_\ell \left( 1 -
        {\text{e}}^{-\tau_{\ell,\text{path}}} \right), \f] between the scattered fraction
        \f$L_\ell^{\text{sca}}\f$ and the \f$N\f$ absorbed fractions \f$L_{\ell,n}^{\text{abs}}\f$.
        When there is only one single dust component, the optical properties of the dust are the
        same everywhere along the path, and we easily have \f[ L_\ell^{\text{sca}} = \varpi_\ell\,
        L_\ell \left( 1 - {\text{e}}^{-\tau_{\ell,\text{path}}} \right) \f] with \f[ \varpi_\ell =
        \frac{ \kappa_\ell^{\text{sca}} }{ \kappa_\ell^{\text{ext}} } \f] the scattering albedo of
        the dust, and \f[ L_{\ell,n}^{\text{abs}} = (1-\varpi_\ell)\, L_\ell \left(
        {\text{e}}^{-\tau_{\ell,n-1}} - {\text{e}}^{-\tau_{\ell,n}} \right), \f] with
        \f$\tau_{\ell,n}\f$ the optical depth measured from the initial position of the path until
        the exit point of the \f$n\f$'th dust cell along the path (this quantity is also stored in
        the PhotonPacket object). It is straightforward to check that the sum of escaped,
        scattered and absorbed luminosities is \f[ L_\ell^{\text{esc}} + L_\ell^{\text{sca}} +
        \sum_{n=0}^{N-1} L_{\ell,n}^{\text{abs}} = L_\ell, \f] as desired. When there is more than
        one dust component, the dust properties are no longer uniform in every cell and things
        become a bit more messy. The expression for \f$L_{\ell,n}^{\text{abs}}\f$ can be readily
        expanded from the corresponding expression in the case of a single dust component, \f[
        L_{\ell,n}^{\text{abs}} = (1-\varpi_{\ell,n})\, L_\ell \left( {\text{e}}^{-\tau_{\ell,n-1}}
        - {\text{e}}^{-\tau_{\ell,n}} \right). \f] The only change is that the global albedo
        \f$\varpi_\ell\f$ is now replaced by a local albedo \f$\varpi_{\ell,n}\f$. It is found as
        the weighted mean of the albedo of the different dust components in the \f$n\f$'th dust
        cell along the path, \f[ \varpi_{\ell,n} = \frac{ \sum_h \kappa_{\ell,h}^{\text{sca}}\,
        \rho_{h,n} }{ \sum_h \kappa_{\ell,h}^{\text{ext}}\, \rho_{h,n} } \f] with \f$\rho_{h,n}\f$
        the density of the \f$h\f$'th dust component in the \f$n\f$'th cell. To calculate
        \f$L_\ell^{\text{sca}}\f$, we also have to take into account that the albedo varies from
        cell to cell, \f[ L_\ell^{\text{sca}} = L_\ell\, \sum_{n=0}^{N-1} \varpi_{\ell,n} \left(
        {\text{e}}^{-\tau_{\ell,n-1}} - {\text{e}}^{-\tau_{\ell,n}} \right). \f] Also in this case,
        it is easy to see that \f[ L_\ell^{\text{esc}} + L_\ell^{\text{sca}} + \sum_{n=0}^{N-1}
        L_{\ell,n}^{\text{abs}} = L_\ell. \f] */
    void simulateEscapeAndAbsorption(PhotonPacket* pp, bool storeAbsorption);

    /** This function determines the next scattering location of a photon packet and the simulates
        the propagation to this position. Given the total optical depth along the path of the
        photon packet \f$\tau_{\ell,\text{path}}\f$ (this quantity is stored in the PhotonPacket
        object provided as an input parameter of this function), the appropriate probability
        distribution for the covered optical depth is an exponential probability
        distribution cut off at \f$\tau_{\ell,\text{path}}\f$. Properly normalized, it reads as
        \f[ p(\tau_\ell) = \frac{{\text{e}}^{-\tau_\ell}}
        {1-{\text{e}}^{-\tau_{\ell,\text{path}}}} \f] where the range of \f$\tau_\ell\f$ is limited
        to the interval \f$[0,\tau_{\ell,\text{path}}]\f$. Instead of generating a
        random optical depth \f$\tau_\ell\f$ directly from this distribution, we use the biasing
        technique in order to cover the entire allowed optical depth range \f$[0,\tau_{\ell,
        \text{path}}]\f$ more uniformly. As the biased probability distribution, we use a linear
        combination between an exponential distribution and a uniform distribution, with a parameter
        \f$\xi\f$ setting the relative importance of the uniform part. In formula form, \f[
        q(\tau_\ell) = (1-\xi)\, \frac{ {\text{e}}^{-\tau_\ell} }
        { 1-{\text{e}}^{-\tau_{\ell,\text{path}}} } + \frac{\xi}{\tau_{\ell,\text{path}}}. \f] A
        random optical depth from this distribution is readily determined. Since we use biasing, the
        weight, or correspondingly the luminosity, of the photon packet needs to be adjusted with
        a bias factor \f$p(\tau_\ell)/q(\tau_\ell)\f$. Finally, the randomly determined optical
        depth is converted to a physical path length \f$s\f$, and the photon packet is
        propagated over this distance. */
    void simulatePropagation(PhotonPacket* pp);

    /** This function simulates the peel-off of a photon packet before a scattering event. This
        means that, just before a scattering event, we create peel-off or shadow photon packets,
        one for every instrument in the instrument system, that we force to propagate in the
        direction of the observer(s) instead of in the propagation direction \f${\bf{k}}\f$
        determined randomly by the scattering process. Each peel-off photon packet has the same
        characteristics as the original photon packet, apart from three differences. The first one
        is, obviously, that the propagation direction is altered to the direction
        \f${\bf{k}}_{\text{obs}}\f$ towards the observer. The second difference is that we have to
        alter the luminosity of the photon packet to compensate for this change in propagation
        direction. This compensation is necessary because the scattering process is anisotropic.
        Since we force the peel-off photon packet to be scattered from the direction
        \f${\bf{k}}\f$ into the direction \f${\bf{k}}_{\text{obs}}\f$, we give it as additional
        weight factor the probability that a photon packet would be scattered into the direction
        \f${\bf{k}}_{\text{obs}}\f$ if its original propagation direction was \f${\bf{k}}\f$. If
        there is only one dust component in the dust system, this weight factor is equal to the
        scattering phase function \f$w= \Phi_\ell({\bf{k}},{\bf{k}}_{\text{obs}})\f$ of the dust
        mixture at the wavelength index of the photon packet. If there are different dust
        components, each with their own density and dust mixture, the appropriate weight factor is
        a weighted mean of the scattering phase functions, \f[ w = \frac{ \sum_h
        \kappa_{\ell,h}^{\text{sca}}\, \rho_{m,h}\, \Phi_{\ell,h}({\bf{k}}, {\bf{k}}_{\text{obs}})
        }{ \sum_h \kappa_{\ell,h}^{\text{sca}}\, \rho_{m,h} }, \f] where \f$\rho_{m,h}\f$ is the
        density of the dust corresponding to the \f$h\f$'th dust component in the dust cell where
        the scattering event takes place, and \f$\kappa_{\ell,h}^{\text{sca}}\f$ and
        \f$\Phi_{\ell,h}\f$ are the scattering coefficient and phase function corresponding to the
        \f$h\f$'th dust component respectively (both evaluated at the wavelength index \f$\ell\f$
        of the photon packet). The third difference is that the polarization state of the peel off
        photon packet is adjusted. If there are multiple dust components, the weight factors
        described above are used not just for the luminosity but also for the components of the
        Stokes vector. For each instrument in the instrument system, the function creates such a
        peel-off photon packet and feeds it to the instrument. The first argument specifies the
        photon packet that was just emitted; the second argument provides a placeholder peel off
        photon packet for use by the function. */
    void peelOffScattering(const PhotonPacket* pp, PhotonPacket* ppp);

    /** This function simulates a scattering event of a photon packet. Most of the properties of
        the photon packet remain unaltered, including the position and the luminosity. The
        properties that change are the number of scattering events experienced by the photon
        packet (this is obviously increased by one) the propagation direction, which is generated
        randomly, and the polarization state. If there is only one dust component, the propagation
        direction and polarization state are obtained from the scattering phase function. If there
        are several components, the function first generates a random dust component (where the
        relative weight of each dust component is equal to \f[ w_h = \frac{
        \kappa_{\ell,h}^{\text{sca}}\, \rho_{h,m} }{ \sum_{h'} \kappa_{\ell,h'}^{\text{sca}}\,
        \rho_{h',m} }, \f] where \f$\rho_{m,h}\f$ is the density of the dust corresponding to the
        \f$h\f$'th dust component in the dust cell where the scattering event takes place, and
        \f$\kappa_{\ell,h}^{\text{sca}}\f$ is the scattering coefficient corresponding to the
        \f$h\f$'th dust component respectively. When a random dust component is generated, a random
        propagation direction and polarization state are obtained from the corresponding scattering
        phase function. */
    void simulateScattering(PhotonPacket* pp);

    //======================== Data Members ========================

private:
    // non-discoverable simulation items
    Configuration* _config{ new Configuration(this) };

    // data members used by the XXXprogress() functions in this class
    string _segment;               // a string identifying the photon shooting segment for use in the log message
    size_t _numTotal;              // the total number of photon packets to be processed for this segment
};

////////////////////////////////////////////////////////////////////

#endif
