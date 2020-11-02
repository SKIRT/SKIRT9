/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MONTECARLOSIMULATION_HPP
#define MONTECARLOSIMULATION_HPP

#include "Configuration.hpp"
#include "Cosmology.hpp"
#include "InstrumentSystem.hpp"
#include "MediumSystem.hpp"
#include "ProbeSystem.hpp"
#include "Simulation.hpp"
#include "SourceSystem.hpp"
class SecondarySourceSystem;

//////////////////////////////////////////////////////////////////////

/** The MonteCarloSimulation class is the top-level class describing a SKIRT simulation. Running a
    Monte Carlo simulation with SKIRT essentially comes down to constructing an instance of the
    MonteCarloSimulation class and invoking the setupAndRun() function on it. The
    MonteCarloSimulation class holds the source, media, instrument and probe systems, implements
    the core aspects of the photon packet life-cycle, and manages the iterative processes in the
    simulation (including phases, iterations and segments).

    The user-configurable \em SimulationMode enumeration sets the overall simulation mode. It
    determines the simulation's wavelength regime (oligochromatic or panchromatic) and the required
    simulation segments (e.g., extinction only or including secondary emission). The
    runSimulation() function performs all segments of the simulation, calling the
    runPrimaryEmission(), runDustSelfAbsorptionPhase(), and runSecondaryEmission() functions as
    required. These functions in turn invoke the performLifeCycle() function to trace photon
    packets through their complete life cycle, including emission, multiple forced scattering
    events, peel-off towards the instruments, and registration of the contribution to the radiation
    field in each spatial cell crossed.

    The MonteCarloSimulation class also holds the non-discoverable \em config property, which is
    automatically set to an instance of the Configuration class. The setup() function of the config
    object is invoked at the very early stages of overall simulation setup, so that it can
    initialize its internal state to reflect the simulation configuration. As a result, it is safe
    for other simulation items to retrieve information from the config object during setup. */
class MonteCarloSimulation : public Simulation
{
    /** The enumeration type indicating the simulation mode, which determines the overall structure
        of the simulation and its capabilities. The choice made for the simulation mode has a
        significant impact on which options are allowed or required in the simulation's
        configuration. For example, some basic simulation modes support just primary sources
        without any media.

        An important aspect determined by the simulation mode is the simulation's wavelength
        regime, which can be oligochromatic or panchromatic. Oligochromatic simulations use just a
        few pre-defined, discrete wavelengths. They do no support kinematics (moving sources and/or
        media) because the wavelengths cannot shift away from the pre-defined values, and they do
        not support secondary emission by the transfer medium because the radiation field must be
        known across a wide spectrum to calculate the medium state and the resulting emission.
        Panchromatic simulations use a continuous range of wavelengths, lifting these limitations.

        For panchromatic simulations, the simulation mode further determines fundamental choices
        such as whether to include secondary emission, and whether to iterate over the radiation
        field state to obtain self-consistent results taking into account, for example, dust
        self-absorption.

        The "extinction-only" simulation modes calculate the extinction of the primary radiation
        through the configured media, including the effects of absorption and scattering. There is
        no secondary emission, so these modes are meaningful only for wavelengths at which
        secondary sources (radiation from the media) can be neglected, i.e. in the ultraviolet,
        optical and near-infrared. Also, the media state is constant, i.e. it is initialized from
        the properties defined in the input model (density distributions, material properties) and
        is never updated. As a result, there is no need to store the radiation field during the
        photon packet life cycle. However, there is user-configurable option to store the radiation
        field anyway so that it can be probed for output.

        The "dust emission" simulation modes (which require a panchromatic wavelength range that
        includes optical to far-infrared regimes) include secondary emission from dust, in addition
        to the effects of absorption and scattering. In these modes, the simulation keeps track of
        the radation field (to calculate the dust emission) and optionally performs iterations to
        self-consistently calculate the effects of dust self-absorption. */
    ENUM_DEF(SimulationMode, OligoNoMedium, OligoExtinctionOnly, NoMedium, ExtinctionOnly, DustEmission,
             DustEmissionWithSelfAbsorption, LyaWithDustExtinction)
        ENUM_VAL(SimulationMode, OligoNoMedium, "No medium - oligochromatic regime (a few discrete wavelengths)")
        ENUM_VAL(SimulationMode, OligoExtinctionOnly,
                 "Dust extinction only - oligochromatic regime (a few discrete wavelengths)")
        ENUM_VAL(SimulationMode, NoMedium, "No medium (primary sources only)")
        ENUM_VAL(SimulationMode, ExtinctionOnly, "Dust extinction only (no secondary emission)")
        ENUM_VAL(SimulationMode, DustEmission, "With secondary emission from dust")
        ENUM_VAL(SimulationMode, DustEmissionWithSelfAbsorption,
                 "With secondary emission from dust and iterations for dust self-absorption")
        ENUM_VAL(SimulationMode, LyaWithDustExtinction, "Lyman-alpha line transfer and dust extinction")
    ENUM_END()

    ITEM_CONCRETE(MonteCarloSimulation, Simulation, "a Monte Carlo simulation")

        PROPERTY_ENUM(simulationMode, SimulationMode, "the overall simulation mode")
        ATTRIBUTE_DEFAULT_VALUE(simulationMode, "ExtinctionOnly")
        ATTRIBUTE_INSERT(simulationMode, "simulationModeOligoNoMedium:Oligochromatic,NoMedium;"
                                         "simulationModeOligoExtinctionOnly:Oligochromatic,ExtinctionOnly;"
                                         "simulationModeNoMedium:Panchromatic,NoMedium;"
                                         "simulationModeExtinctionOnly:Panchromatic,ExtinctionOnly;"
                                         "simulationModeDustEmission:Panchromatic,DustEmission,Emission,RadiationField;"
                                         "simulationModeDustEmissionWithSelfAbsorption:"
                                         "Panchromatic,DustEmission,Emission,RadiationField,DustSelfAbsorption;"
                                         "simulationModeLyaWithDustExtinction:Lya,Panchromatic,ExtinctionOnly")

        PROPERTY_ITEM(cosmology, Cosmology, "the cosmology parameters")
        ATTRIBUTE_DEFAULT_VALUE(cosmology, "LocalUniverseCosmology")
        ATTRIBUTE_DISPLAYED_IF(cosmology, "Level2")

        PROPERTY_DOUBLE(numPackets, "the default number of photon packets launched per simulation segment")
        ATTRIBUTE_MIN_VALUE(numPackets, "[0")
        ATTRIBUTE_MAX_VALUE(numPackets, "1e19]")
        ATTRIBUTE_DEFAULT_VALUE(numPackets, "1e6")

        PROPERTY_ITEM(sourceSystem, SourceSystem, "the source system")
        ATTRIBUTE_DEFAULT_VALUE(sourceSystem, "SourceSystem")

        PROPERTY_ITEM(mediumSystem, MediumSystem, "the medium system")
        ATTRIBUTE_DEFAULT_VALUE(mediumSystem, "!NoMedium:MediumSystem;")
        ATTRIBUTE_REQUIRED_IF(mediumSystem, "!NoMedium")

        PROPERTY_ITEM(instrumentSystem, InstrumentSystem, "the instrument system")
        ATTRIBUTE_DEFAULT_VALUE(instrumentSystem, "InstrumentSystem")

        PROPERTY_ITEM(probeSystem, ProbeSystem, "the probe system")
        ATTRIBUTE_DEFAULT_VALUE(probeSystem, "ProbeSystem")

    ITEM_END()

    /** \fn numPackets
        The number of photon packets is specified as a double-precision floating point number
        rather than as a 64-bit integer to avoid implementing yet another discoverable property
        type. As a side benefit, one can use exponential notation to specify a large number of
        photon packets. Also, note that a double can exactly represent all integers up to 9e15. The
        maximum number of photon packets is somewhat arbitrarily set to 1e19 because that number is
        close to the maximum number representable with a 64-bit unsigned integer. */

    //============= Construction - Setup - Destruction =============

protected:
    /** This function performs setup for the complete simulation hierarchy. It calls the regular
        setup() function and notifies the probe system when setup has been completed. */
    void setupSimulation() override;

    /** This function performs initial setup for the MonteCarloSimulation object. For example,
        it constructs a SecondarySourceSystem object if the simulation configuration requires
        secondary emission. */
    void setupSelfBefore() override;

    //======== Getters for Non-Discoverable Properties =======

public:
    /** Returns the Configuration object for this simulation hierarchy. */
    Configuration* config() const;

    //======================== Other Functions =======================

protected:
    /** This function actually runs the simulation, assuming setup has been completed. It
        consecutively performs all segments of the simulation, including primary source segments,
        secondary source segments, and any required iterations. For each segment, it tells the
        source system to prepare for launching photon packets in serial code, and then it causes
        photon packets to be launched (and traced through their life cycles) by calling the
        performLifeCycle() function in appropriately parallelized code depending on the run-time
        environment and the command-line options. After each of the segments but the last one, the
        function also tells the medium system to synchronize the radiation field between processes
        (in a multi-process environment). */
    void runSimulation() override;

private:
    /** This function runs the primary source emission segment. It implements a parallelized loop
        that iterates over the configured number of photon packets and drives the photon packet
        life cycle for each. The primary emission segment includes peel-off and records radiation
        field contributions if the configuration requires it (e.g., because secondary emisison must
        be calculated, or because the user configured probes to directly output radiation field
        information.) */
    void runPrimaryEmission();

    /** This function implements primary source emission when the simulation includes a dynamic
        medium state. In that case, the primary emission phase self-consistently calculates the
        radiation field taking into account one or more recipes for updating the medium state as a
        function of the radiation field. See the DynamicStateRecipe class for more information.

        Specifically, this function implements a number of iterations, where each iteration tracks
        a segment of photon packets through the medium and subsequently updates the medium state
        based on the newly established radiation field. The number of photon packets launched for
        each iteration can be configured as a multiplication factor on the default number of
        primary photon packets. The minimum and maximum number of iterations can also be specified
        as configuration options. Within these limits, the actual number of iterations performed is
        determined by convergence criteria defined by the configured dynamic medium state
        recipe(s).

        After the dynamic state iterations have completed (with or without convergence), this
        function performs a final segment of primary emission photon launching that now includes
        peel-off towards the instruments. */
    void runPrimaryEmissionWithDynamicState();

    /** This function runs the dust self-absorption phase. This phase includes a series of
        intermediate secondary source emission segments in an iteration to self-consistently
        calculate the radiation field, taking into account the fraction of dust emission absorbed
        by the dust itself. The intermediate secondary emission segments do not perform peel-off
        towards the instruments (because the radiation field is not yet converged). They record
        radiation field contributions in a seperate table, so that the radiation field resulting
        from primary emission remains untouched and can be reused.

        The minimum and maximum number of iterations can be specified as configuration options.
        Within these limits, the actual number of iterations performed is determined by convergence
        criteria which can also be specified as configuration options. Convergence is reached (and
        the function exits) when (a) the absorbed dust luminosity is less than a given fraction of
        the absorbed stellar luminosity, \em OR (b) the absorbed dust luminosity has changed by
        less than a given fraction compared to the previous iteration. */
    void runDustSelfAbsorptionPhase();

    /** This function runs the final secondary source emission segment. It implements a
        parallelized loop that iterates over the configured number of photon packets and drives the
        photon packet life cycle for each. The final secondary emission segment includes peel-off
        but does not record radiation field contributions, because the radiation field is assumed
        to be converged (or simply immutable). */
    void runSecondaryEmission();

    /** In a multi-processing environment, this function logs a message and waits for all processes
        to finish the work (i.e. it places a barrier). The string argument is included in the log
        message to indicate the scope of work that is being finished. If there is only a single
        process, the function does nothing. */
    void wait(string scope);

    /** This function initializes the progress counter used in logprogress() for the specified
        segment and logs the number of photon packets to be processed. */
    void initProgress(string segment, size_t numTotal);

    /** This function logs a progress message for the segment specified in the initprogress()
        function if the previous message was issued at least some time interval ago. The function
        must be called regularly while processing photon packets. The argument specifies the number
        of photon packets processed. */
    void logProgress(size_t numDone);

    /** This function launches the specified chunk of photon packets from primary or secondary
        sources, and it implements the complete life-cycle for each of these photon packets. This
        includes emission and multiple forced scattering events, and, if requested, the
        corresponding peel-off towards the instruments, and registration of the contribution to the
        radiation field in each spatial cell crossed.

        A photon packet is born when it is emitted by either the primary or the secondary source
        system. Immediately after birth, peel-off photon packets are created and launched towards
        the instruments (one for each instrument). If the simulation contains one or more media,
        the photon packet now enters a cycle which consists of different steps. First, all details
        of the path of the photon packet through the medium are calculated and stored. Based on
        this information, the contribution of the photon packet to the radiation field in each
        spatial cell can be stored. The packet is then propagated to a new (random) scattering
        position, and its weight (luminosity) is adjusted for the escape fraction and the absorbed
        energy. If the photon packet is still sufficiently luminous, scattering peel-off photon
        packets are created and launched towards each instrument, and the actual scattering event
        is simulated. Finally, the loop repeats itself. It is terminated only when the photon
        packet has lost a substantial part of its original luminosity (and hence becomes
        irrelevant).

        The first two arguments of this function specify the range of photon packet history indices
        to be handled. The \em primary flag is true to launch from primary sources, false for
        secondary sources. The \em peel flag indicates whether peeloff photon packets should be
        sent towards the instruments. The \em store flag indicates whether the contribution to the
        radiation field should be stored. */
    void performLifeCycle(size_t firstIndex, size_t numIndices, bool primary, bool peel, bool store);

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
        probability that it is emitted in any other direction. If the source has a nonzero velocity,
        the wavelength of the peel-off photon packet is Doppler-shifted for the new direction.

        The first argument specifies the photon packet that was just emitted; the second argument
        provides a placeholder peel off photon packet for use by the function. */
    void peelOffEmission(const PhotonPacket* pp, PhotonPacket* ppp);

    /** This function stores the contribution of the specified photon packet to the radiation field
        in the cells crossed by the packet's path. The function assumes that both the geometric and
        optical depth information for the photon packet's path have been set; if this is not the
        case, the behavior is undefined.

        The radiation field is discretized on the spatial grid of the simulation (held by the
        MediumSystem object) and on a wavelength grid provided specifically for this purpose (see
        the Configuration::radiationFieldWLG() function for more information). Locating the
        appropriate spatial bin is trivial because each segment in the photon packet's path stores
        the index of the cell being crossed. The wavelength bin is derived from the photon packet's
        perceived wavelength in the cell under consideration, taking into account the velocity
        of the medium in that cell.

        For each segment \f$n\f$ in the photon packet's path, this function first determines the
        corresponding spatial/wavelength bin as described above. To the contents of that bin, the
        function adds the product of the mean luminosity carried by the photon packet along the
        segment and the distance covered by the segment, or \f[ (L\Delta s)_n = L_{\text{mean},n}
        \,(\Delta s)_n = L\,\text{lnmean} \left(\text{e}^{-\tau_{n-1}}, \text{e}^{-\tau_n}\right)
        \,(\Delta s)_n \f] where \f$L\f$ is the luminosity carried by the photon packet,
        \f$\text{lnmean}()\f$ is the logarithmic mean, \f$\tau_{n-1}\f$ and \f$\tau_n\f$ represent
        the cumulative optical depth at the start and end of the segment, and \f$(\Delta s)_n\f$ is
        the distance covered by the segment. Using the logarithmic mean assumes an exponential
        behavior of the exinction with distance within the segment.

        Once this information has been accumulated for all photon packets launched during a segment
        of the simulation, the mean intensity of the radiation field in each spatial/wavelength bin
        can be calculated using \f[ (J_\lambda)_{\ell,m} = \frac{ (L\Delta s)_{\ell,m}
        }{4\pi\,V_m\,(\Delta \lambda)_\ell} \f] where \f$\ell\f$ is the index of the wavelength
        bin, \f$(\Delta \lambda)_\ell\f$ is the wavelength bin width, \f$m\f$ is the spatial cell
        index, \f$V_m\f$ is the volume of the cell, and \f$(L\Delta s)_{\ell,m}\f$ has been
        accumulated over all photon packets contributing to the bin. The resulting mean intensity
        \f$J_\lambda\f$ is expressed as an amount of energy per unit of time, per unit of area, per
        unit of wavelength, and per unit of solid angle. */
    void storeRadiationField(const PhotonPacket* pp);

    /** This function determines the next scattering location of a photon packet in a photon life
        cycle with forced scattering and simulates its propagation to that position. The function
        assumes that both the geometric and optical depth information for the photon packet's path
        have been set; if this is not the case, the behavior is undefined. This function proceeds
        in a number of steps as outlined below.

        <b>Total optical depth</b>

        We first determine the total optical depth \f$\tau_\text{path}\f$ of the photon packet's
        path. Because the path has been calculated until the edge of the simulation's spatial grid,
        \f$\tau_\text{path}\f$ is equal to the cumulative optical depth at the end of the last
        segment in the path.

        <b>%Random optical depth</b>

        We then randomly generate the optical depth \f$\tau\f$ at the interaction site from an
        appropriate distribution. Given the total optical depth along the path of the photon packet
        \f$\tau_\text{path}\f$, the appropriate probability distribution for the covered optical
        depth is an exponential probability distribution cut off at \f$\tau_\text{path}\f$.
        Properly normalized, it reads as \f[ p(\tau) = \frac{{\text{e}}^{-\tau}}
        {1-{\text{e}}^{-\tau_\text{path}}} \f] where the range of \f$\tau\f$ is limited to the
        interval \f$[0,\tau_\text{path}]\f$. Instead of generating a random optical depth
        \f$\tau\f$ directly from this distribution, we use the biasing technique in order to cover
        the entire allowed optical depth range \f$[0,\tau_\text{path}]\f$ more uniformly. As the
        biased probability distribution, we use a linear combination between an exponential
        distribution and a uniform distribution, with a parameter \f$\xi\f$ setting the relative
        importance of the uniform part. In formula form, \f[ q(\tau) = (1-\xi)\, \frac{
        {\text{e}}^{-\tau} } { 1-{\text{e}}^{-\tau_\text{path}} } + \frac{\xi}{\tau_\text{path}}.
        \f] A random optical depth from this distribution is readily determined. Since we use
        biasing, the weight, or correspondingly the luminosity, of the photon packet needs to be
        adjusted with a bias factor \f$p(\tau)/q(\tau)\f$.

        <b>Interaction point</b>

        Now that the optical depth \f$\tau\f$ at the interaction site has been (randomly) chosen,
        we determine the physical position of the interaction point along the path. This is
        accomplished in two steps: a binary search among the path segments to determine the segment
        (or cell) "containing" the given cumulative optical depth, and subsequent linear
        interpolation within the cell assuming exponential behavior of the extinction.

        <b>Albedo</b>

        We calculate the scattering albedo \f$\varpi\f$ of the medium at the interaction point (or
        more precisely, for the spatial cell containing the interaction point).

        <b>Weight adjustment</b>

        We adjust the weight of the photon packet to compensate for the escaped and absorbed
        portions of the luminosity. More precisely, the weight is multiplied by the scattered
        fraction, i.e. the fraction of the luminosity that does not escape and does not get
        absorbed, \f[ f_\text{sca} = (1-f_\text{esc}) \,\varpi = (1-\text{e}^{-\tau_\text{path}})
        \,\varpi. \f]

        <b>Advance position</b>

        Finally we advance the initial position of the photon packet to the interaction point. This
        last step invalidates the photon packet's path (including geometric and optical depth
        information). The packet is now ready to be scattered into a new direction. */
    void simulateForcedPropagation(PhotonPacket* pp);

    /** This function simulates the propagation of a photon packet to the next scattering location
        in a photon life cycle without forced scattering. The function assumes that the next
        scattering location for the photon packet's path has already been set; if this is not the
        case, the behavior is undefined. This function proceeds in a number of steps as outlined
        below.

        <b>%Random optical depth</b>

        Because the photon packet is allowed to escape the model, we randomly generate the optical
        depth \f$\tau\f$ at the interaction site from the regular exponential distribution. The
        current implementation does not support path length biasing.

        <b>Interaction point</b>

        Now that the optical depth \f$\tau\f$ at the interaction site has been (randomly) chosen,
        we determine the physical position of the interaction point along the path. This is
        accomplished by generating the path segments for crossed cells one by one and calculating
        the corresponding cumulative optical depth on the fly. Once the cumulative optical depth at
        the exit point of a segment exceeds the interaction optical depth, we know that the
        interaction point must be within this segment. The physical interaction point is then
        obtained through linear interpolation within the segment, assuming exponential behavior of
        the extinction. If the cumulative optical depth of the path never exceeds the interaction
        optical depth, the photon packet escapes and is terminated.

        <b>Albedo</b>

        We calculate the scattering albedo \f$\varpi\f$ of the medium at the interaction point (or
        more precisely, for the spatial cell containing the interaction point).

        <b>Weight adjustment</b>

        We adjust the weight of the photon packet to compensate for the absorbed portion of the
        luminosity. More precisely, the weight is multiplied by the scattered fraction, i.e. \f$
        f_\text{sca} = \varpi \f$.

        <b>Advance position</b>

        Finally we advance the initial position of the photon packet to the interaction point. This
        last step invalidates the photon packet's path (including geometric and optical depth
        information). The packet is now ready to be scattered into a new direction. */
    void simulateNonForcedPropagation(PhotonPacket* pp);

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
        account the bulk velocity of the medium.

        The first argument to this function specifies the photon packet that is about to be
        scattered; the second argument provides a placeholder peel off photon packet for use by the
        function. */
    void peelOffScattering(PhotonPacket* pp, PhotonPacket* ppp);

    //======================== Data Members ========================

private:
    // non-discoverable simulation items
    Configuration* _config{new Configuration(this)};
    SecondarySourceSystem* _secondarySourceSystem{nullptr};  // constructed only when there is secondary emission

    // data members used by the XXXprogress() functions in this class
    string _segment;  // a string identifying the photon shooting segment for use in the log message
};

////////////////////////////////////////////////////////////////////

#endif
