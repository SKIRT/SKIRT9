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

    /** This function launches the specified chunk of photon packets from primary sources. It
        implements the complete photon packet life-cycle, including emission and multiple forced
        scattering events, as well as the corresponding peel-off photon packets towards the
        instruments. If required for the configuration of the simulation, the function also
        registers the contribution of the photon packet to the radiation field in each spatial
        cell crossed by its path. */
    void doPrimaryEmissionChunk(size_t firstIndex, size_t numIndices);

    /** This function implements the peel-off of a photon packet after an emission event. This
        means that we create a peel-off photon packet for every instrument in the instrument
        system, which is forced to propagate in the direction of the observer instead of in the
        propagation direction determined randomly by the emission process. Each peel-off photon
        packet is subsequently fed into its target instrument for detection.

        A peel-off photon packet has the same characteristics as the original photon packet, except
        that the propagation direction is altered from the emission direction \f${\bf{k}}\f$ to the
        direction of the observer \f${\bf{k}}_{\text{obs}}\f$. For anisotropic emission, a weight
        factor is applied to the luminosity to compensate for the fact that the probability that a
        photon packet would have been emitted towards the observer is not the same as the
        probability that it is emitted in any other direction. If the source has a bulk velocity,
        the wavelength of the peel-off photon packet is Doppler-shifted for the new direction.

        The first argument specifies the photon packet that was just emitted; the second argument
        provides a placeholder peel off photon packet for use by the function. */
    void peelOffEmission(const PhotonPacket* pp, PhotonPacket* ppp);

    /** This function simulates the escape from the system and the absorption by the transfer
        medium of a fraction of the luminosity of a photon packet. It actually splits the
        luminosity \f$L\f$ of the photon packet in \f$N+2\f$ different parts, with \f$N\f$ the
        number of spatial cells along its path through the medium: a part \f$L^{\text{esc}}\f$ that
        escapes from the system, a part \f$L^{\text{sca}}\f$ that is scattered before it leaves the
        system, and \f$N\f$ parts \f$L_n^{\text{abs}}\f$, each of them corresponding to the
        fraction that is absorbed in a cell along the path. The part that scatters is the part of
        the photon packet that survives and continues in the photon packet's life cycle. So we
        alter the luminosity of the photon packet from \f$L\f$ to \f$L^{\text{sca}}\f$.

        In the analysis below, we drop the wavelength-dependency of the material properties from
        the notation. Given the total optical depth along the path of the photon packet
        \f$\tau_\text{path}\f$ (this quantity is stored in the PhotonPacket object provided as an
        input parameter to this function), the fraction of the luminosity that escapes from the
        system without any interaction is \f[ L^{\text{esc}} = L\, {\text{e}}^{-\tau_\text{path}}.
        \f] Remains to subdivide the remainder of the initial luminosity, \f[ L \left( 1 -
        {\text{e}}^{-\tau_\text{path}} \right), \f] between the scattered fraction
        \f$L^{\text{sca}}\f$ and the \f$N\f$ absorbed fractions \f$L_n^{\text{abs}}\f$. Taking into
        account that the medium properties may differ between spatial cells, \f$L_n^{\text{abs}}\f$
        can written as \f[ L_n^{\text{abs}} = (1-\varpi_n)\, L \left( {\text{e}}^{-\tau_{n-1}} -
        {\text{e}}^{-\tau_n} \right). \f] The local albedo \f$\varpi_n\f$ is found as the weighted
        mean of the albedo of the different medium components in the \f$n\f$'th cell along the
        path, \f[ \varpi_n = \frac{ \sum_h \varsigma_h^{\text{sca}}\, n_{h,n} }{ \sum_h
        \varsigma_h^{\text{ext}}\, n_{h,n} } \f] with \f$n_{h,n}\f$ the number density of the
        \f$h\f$'th medium component in the \f$n\f$'th cell. In a similar fashion, we can write
        \f$L^{\text{sca}}\f$ as \f[ L^{\text{sca}} = L\, \sum_{n=0}^{N-1} \varpi_n \left(
        {\text{e}}^{-\tau_{n-1}} - {\text{e}}^{-\tau_n} \right). \f] It is easy to see that \f[
        L^{\text{esc}} + L^{\text{sca}} + \sum_{n=0}^{N-1} L_n^{\text{abs}} = L. \f]

        The albedo \f$\varpi_n\f$ in the above analysis is calculated for each spatial cell at the
        wavelength perceived by the medium in that cell taking into account the local bulk
        velocity.

        The first argument specifies the photon packet that needs to be handled; its luminosity is
        updated as described above. If the second argument is set to true, the function also
        registers the contribution of the photon packet to the radiation field in each spatial cell
        crossed by its path. */
    void simulateEscapeAndAbsorption(PhotonPacket* pp, bool storeAbsorption);

    /** This function determines the next scattering location of a photon packet and the simulates
        the propagation to this position. In the analysis below, we drop the wavelength-dependency
        of the optical depth from the notation. Given the total optical depth along the path of the
        photon packet \f$\tau_\text{path}\f$ (this quantity is stored in the PhotonPacket object
        provided as an input parameter to this function), the appropriate probability distribution
        for the covered optical depth is an exponential probability distribution cut off at
        \f$\tau_\text{path}\f$. Properly normalized, it reads as \f[ p(\tau) =
        \frac{{\text{e}}^{-\tau}} {1-{\text{e}}^{-\tau_\text{path}}} \f] where the range of
        \f$\tau\f$ is limited to the interval \f$[0,\tau_\text{path}]\f$. Instead of generating a
        random optical depth \f$\tau\f$ directly from this distribution, we use the biasing
        technique in order to cover the entire allowed optical depth range
        \f$[0,\tau_\text{path}]\f$ more uniformly. As the biased probability distribution, we use a
        linear combination between an exponential distribution and a uniform distribution, with a
        parameter \f$\xi\f$ setting the relative importance of the uniform part. In formula form,
        \f[ q(\tau) = (1-\xi)\, \frac{ {\text{e}}^{-\tau} } { 1-{\text{e}}^{-\tau_\text{path}} } +
        \frac{\xi}{\tau_\text{path}}. \f] A random optical depth from this distribution is readily
        determined. Since we use biasing, the weight, or correspondingly the luminosity, of the
        photon packet needs to be adjusted with a bias factor \f$p(\tau)/q(\tau)\f$. Finally, the
        randomly determined optical depth is converted to a physical path length \f$s\f$, and the
        photon packet is propagated over this distance. */
    void simulatePropagation(PhotonPacket* pp);

    /** This function simulates the peel-off of a photon packet before a scattering event. This
        means that, just before a scattering event, we create a peel-off photon packet for every
        instrument in the instrument system, which is forced to propagate in the direction of the
        observer instead of in the propagation direction determined randomly by the scattering
        process. Each peel-off photon packet is subsequently fed into its target instrument for
        detection.

        A peel-off photon packet has the same characteristics as the original photon packet, apart
        from four differences. The first one is, obviously, that the propagation direction is
        altered to the direction towards the observer. The second difference is that we have to
        alter the luminosity of the photon packet to compensate for this change in propagation
        direction, because the scattering process is anisotropic. The third difference is that the
        polarization state of the peel-off photon packet is adjusted. And the last difference is
        that the wavelength of the peel-off photon packet is properly Doppler-shifted taking into
        account the bulk velocity of the medium. We discuss some of these changes in more detail
        below. In this analysis, we drop the wavelength-dependency of the material properties from
        the notation.

        Since we force the peel-off photon packet to be scattered from the direction \f${\bf{k}}\f$
        into the direction \f${\bf{k}}_{\text{obs}}\f$, the corresponding biasing weight factor is
        given by the probability that a photon packet would be scattered into the direction
        \f${\bf{k}}_{\text{obs}}\f$ if its original propagation direction was \f${\bf{k}}\f$. If
        there is only one medium component in the system, this weight factor is equal to the value
        of the scattering phase function \f$w= \Phi({\bf{k}},{\bf{k}}_{\text{obs}})\f$ for that
        medium component. If there are multiple medium components, the weight factor is the
        weighted mean of the scattering phase function values, \f[ w = \frac{ \sum_h
        \varsigma_h^{\text{sca}}\, n_{m,h}\, \Phi_h({\bf{k}}, {\bf{k}}_{\text{obs}}) }{ \sum_h
        \varsigma_h^{\text{sca}}\, n_{m,h} }, \f] where \f$n_{m,h}\f$ is the number density of the
        medium corresponding to the \f$h\f$'th component in the cell where the scattering event
        takes place, and \f$\varsigma_h^{\text{sca}}\f$ and \f$\Phi_h\f$ are the scattering cross
        section and phase function corresponding to the \f$h\f$'th component respectively.

        Evaluation of the phase function depends on the scattering mode supported by each medium's
        material mix. For the most basic mode, the material mix provides a value for the scattering
        asymmetry parameter \f$g=\left<\cos\theta\right>\f$. A value of \f$g=0\f$ corresponds to
        isotropic scattering. Other values \f$-1\le g\le 1\f$ are substituted in the
        Henyey-Greenstein phase function, \f[ \Phi(\cos\theta) = \frac{1-g^2}
        {(1+g^2-2g\cos\theta)^{3/2}}. \f] For other scattering modes, the phase function provided
        by the material mix is invoked instead.

        In case polarization is supported in the current simulation configuration, the polarization
        state of the peel off photon packet is adjusted as well. Note that all media must either
        support polarization or not support it, mixing these support levels is not allowed.
        Compliance with this requirement is verified during setup of the simulation. The adjusted
        Stokes vector for a particular medium component is obtained as follows. The function
        rotates the Stokes vector from the reference direction in the previous scattering plane
        into the peel-off scattering plane, applies the Mueller matrix on the Stokes vector, and
        further rotates the Stokes vector from the reference direction in the peel-off scattering
        plane to the x-axis of the instrument to which the peel-off photon package is headed. If
        there are multiple medium components (all supporting polarization), the weight factors
        described above are used not just for the luminosity but also for the components of the
        Stokes vector.

        The first argument to this function specifies the photon packet that is about to be
        scattered; the second argument provides a placeholder peel off photon packet for use by the
        function. */
    void peelOffScattering(const PhotonPacket* pp, PhotonPacket* ppp);

    /** This function simulates a scattering event of a photon packet. Most of the properties of
        the photon packet remain unaltered, including the position and the luminosity. The
        properties that change are the number of scattering events experienced by the photon packet
        (this is obviously increased by one), the propagation direction, which is generated
        randomly, the polarization state, and the wavelength, which is properly Doppler-shifted
        taking into account the bulk velocity of the medium. In the analysis below, we drop the
        wavelength-dependency of the material properties from the notation.

        If there is only one medium component, the scattering event is governed by the
        corresponding material mix. If there are several components, the function first randomly
        selects a medium component from the list, where the relative weight of each component is
        equal to \f[ w_h = \frac{ \varsigma_h^{\text{sca}}\, n_{h,m} }{ \sum_{h'}
        \varsigma_{\ell,h'}^{\text{sca}}\, n_{h',m} }, \f] where \f$n_{m,h}\f$ is the number
        density of the medium corresponding to the \f$h\f$'th component in the cell where the
        scattering event takes place, and \f$\varsigma_h^{\text{sca}}\f$ is the scattering cross
        section corresponding to the \f$h\f$'th component respectively.

        The remainder of the operation depends on the scattering mode supported by the selected
        medium component's material mix. For the most basic mode, the material mix provides a value
        for the scattering asymmetry parameter \f$g=\left<\cos\theta\right>\f$. For the value
        \f$g=0\f$, corresponding to isotropic scattering, a new direction is generated uniformly on
        the unit sphere. For other values \f$-1\le g\le 1\f$, a scattering angle \f$\theta\f$ is
        sampled from the Henyey-Greenstein phase function, \f[ \Phi(\cos\theta) =
        \frac{1-g^2}{(1+g^2-2g\cos\theta)^{3/2}}. \f] This can be accomplished as follows.
        Substituting \f$\mu=\cos\theta\f$, the probability distribution for \f$\mu\f$ (normalized
        to unity) becomes \f[ p(\mu)\,\text{d}\mu = \frac{1}{2} \,
        \frac{1-g^2}{(1+g^2-2g\mu)^{3/2}} \,\text{d}\mu \qquad -1\leq\mu\leq1 \f] We can use the
        transformation method to sample from this distribution. Given a uniform deviate
        \f$\mathcal{X}\f$, we need to solve \f[ {\mathcal{X}} = \int_{-1}^\mu p(\mu')\,\text{d}\mu'
        \f] Performing the integration and solving for \f$\mu\f$ yields \f[ \cos\theta = \mu =
        \frac{1+g^2-f^2}{2g} \quad\text{with}\quad f=\frac{1-g^2}{1-g+2g {\mathcal{X}}}
        \qquad\text{for}\; g\neq 0 \f] For other scattering modes, a function provided by the
        material mix is invoked instead to obtain a random scattering direction for the photon
        packet.

        In case polarization is supported in the current simulation configuration, the polarization
        state of the photon packet is adjusted as well. Note that all media must either support
        polarization or not support it, mixing these support levels is not allowed. Compliance with
        this requirement is verified during setup of the simulation. The adjusted Stokes vector is
        obtained as follows, again using the randomly selected medium component. After obtaining
        the sampled scattering angles \f$\theta\f$ and \f$\phi\f$ from the material mix, the Stokes
        vector of the photon packet is rotated into the scattering plane and transformed by
        applying the Mueller matrix. Finally, the new direction is computed from the previously
        sampled \f$\theta\f$ and \f$\phi\f$ angles. */
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
