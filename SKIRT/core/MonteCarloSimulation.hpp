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

/** The MonteCarloSimulation class is the top-level class describing a SKIRT simulation. It holds
    the source, media, instrument and probe systems, implements the core aspects of the photon
    packet life-cycle, and manages the iterative processes in the simulation (including phases,
    iterations and segments). Running a Monte Carlo simulation with SKIRT essentially comes down to
    constructing an instance of the MonteCarloSimulation class and invoking the setupAndRun()
    function on it. The runSimulation() function performs all segments of the simulation, calling
    functions such as runPrimaryEmission() and runSecondaryEmission() as required. These functions
    in turn invoke the performLifeCycle() function to trace photon packets through their complete
    life cycle, including emission, scattering events, peel-off towards the instruments, and
    registration of the contribution to the radiation field in each spatial cell crossed.

    The MonteCarloSimulation class also holds the non-discoverable \em config property, which is
    automatically set to an instance of the Configuration class. The setup() function of the config
    object is invoked at the very early stages of overall simulation setup, so that it can
    initialize its internal state to reflect the simulation configuration. As a result, it is safe
    for other simulation items to retrieve information from the config object during setup.

    <b>%Simulation phases, iterations and segments</b>

    In a first simulation phase, SKIRT determines the radiation field (RF) based on just the
    primary emission. In a second phase, it performs secondary emission using spectra calculated
    based on this RF.

    Some SKIRT media allow the density and/or the optical properties of the medium to depend on the
    local RF. The changes to the medium may sufficiently affect the RF to cause a significant
    change in the medium properties itself, so that the RF must be calculated self-consistently by
    iterating over primary and/or secondary emission. We call this a \em dynamic \em medium \em
    state (DMS).

    As a special case, for some SKIRT media the properties depending on the RF are significant only
    for determining the secondary emission spectrum and have no effect during primary emission. As
    a result, there is no need for iteration over primary emission. We call this a \em semi-dynamic
    \em medium \em state (SDMS).

    Secondary emission may sufficiently affect the RF to cause a significant change in the
    secondary emission spectrum itself. If this is the case, the RF must be calculated
    self-consistently by iterating over the secondary emission. We call this \em dynamic \em
    secondary \em emission (DSE).

    A simulation consecutively performs one or more \em segments as part of the relevant \em
    iterations and \em phases. A segment processes a set of photon packets, including emission from
    the relevant sources or media, absorption and scattering by the media, and updating the RF, the
    medium state, and/or the instruments as appropriate. We define the following segment types:

    - \f$\mathbf{P}\f$ -- primary emission: photon packets are launched from the configured
    (primary) source components.

    - \f$\mathbf{S}\f$ -- secondary emission: photon packets are launched from the configured
    medium components that support emission; this includes calculation of the secondary emission
    spectrum for each component based on the RF and the medium state in each spatial cell at the
    start of the segment.

    These segment identifiers are decorated by one or more of the following modifiers:

    - \f$p\f$: peel off a photon packet to the instruments for each emission and scattering event.

    - \f$r\f$: register the RF while tracing photon packets through the medium.

    - \f$(r)\f$: optionally register the RF if requested for probing (not needed for instruments).

    - \f$s\f$: update the DMS based on the RF after processing all photon packets for the segment.

    - \f$(s)\f$: update the SDMS, if any, based on the RF after processing all photon packets for
    the segment.

    <b>%Simulation mode configuration</b>

    The user-configurable \em simulationMode enumeration sets the overall simulation mode. It
    determines the presence or absence of a transfer medium, the wavelength regime (oligochromatic
    or panchromatic), whether there is a secondary emission phase, and if so, which media types
    (dust and/or gas) are emitting photon packets. The Boolean user options \em iterateMediumState
    and \em iterateSecondaryEmission further determine whether the simulation iterates over the
    radiation field to self-consistently calculate a DMS or DSE. The presence of an SDMS is derived
    at run-time from the capabilities advertised by each of the configured medium components.

    An important aspect determined by the simulation mode is the simulation's wavelength regime,
    which can be oligochromatic or panchromatic. Oligochromatic simulations use just a few
    pre-defined, discrete wavelengths. They do no support kinematics (moving sources and/or media)
    because the wavelengths cannot shift away from the pre-defined values, and they do not support
    secondary emission by the transfer medium because the radiation field must be known across a
    wide spectrum to calculate the medium state and the resulting emission. Panchromatic
    simulations use a continuous range of wavelengths, lifting these limitations.

    The "extinction-only" simulation modes calculate the extinction of the primary radiation
    through the configured media, including the effects of absorption and scattering. There is no
    secondary emission, so these modes are meaningful only for wavelengths at which secondary
    sources (radiation from the media) can be neglected. The "emission" simulation modes (which
    require a panchromatic wavelength range) include secondary emission from dust and/or gas, in
    addition to the effects of absorption and scattering.

    The prefix \c Dust-, \c Gas-, or \c DustAndGas- in the emission simulation modes indicates
    which media types are emitting photon packets during secondary emission. Other media types can
    still be configured in the simulation and will contribute extinction during both primary and
    secondary emission.

    The simulation mode and iteration flags are set at the start of the configuration process.
    Their values have a significant impact on which options are allowed or required in the
    simulation's configuration.

    <b>Execution flow</b>

    The following table shows the flow of execution for each simulation mode. \c DMS1 and \c DMS2
    indicate a dynamic medium state during primary and secondary emission, respectively. The arrows
    in the flow diagrams indicate forward progression (\f$\rightarrow\f$) and iteration
    (\f$\longleftarrow\f$).

    | | %Simulation mode | Flow diagram |
    |-|-----------------|--------------|
    | 1 | \c OligoNoMedium | \f$\mathbf{P}^p\f$ |
    | 2 | \c OligoExtinctionOnly | \f$\mathbf{P}^p_{(r)}\f$ |
    | 3 | \c LyaExtinctionOnly | \f$\mathbf{P}^p_{(r)}\f$ |
    | 4 | \c NoMedium | \f$\mathbf{P}^p\f$ |
    | 5 | \c ExtinctionOnly | \f$\mathbf{P}^p_{(r)}\f$ |
    | 6 | ... + \c DMS1 | \f$\overleftarrow{\mathbf{P}_{rs}} \;\rightarrow\; \mathbf{P}^p_{(r)}\f$ |
    | 7 | \c Dust-, \c Gas-, \c DustAndGasEmission | \f$\mathbf{P}^p_{r(s)} \;\rightarrow\; \mathbf{S}^p_{(r)}\f$ |
    | 8 | ... + \c DMS1 | \f$\overleftarrow{\mathbf{P}_{rs}} \;\rightarrow\; \mathbf{P}^p_{r(s)} \;\rightarrow\; \mathbf{S}^p_{(r)}\f$ |
    | 9 | ... + \c DSE | \f$\mathbf{P}^p_{r(s)} \;\rightarrow\; \overleftarrow{\mathbf{S}_{r(s)}} \;\rightarrow\; \mathbf{S}^p_{(r)}\f$ |
    | 10 | ... + \c DSE + \c DMS1 | \f$\overleftarrow{\mathbf{P}_{rs}} \;\rightarrow\; \textbf{}\mathbf{P}^p_{r(s)} \;\rightarrow\; \overleftarrow{\mathbf{S}_{r(s)}} \;\rightarrow\; \mathbf{S}^p_{(r)}\f$ |
    | 11 | ... + \c DSE + \c DMS1 + \c DMS2 | \f$\overleftarrow{\mathbf{P}_{rs}} \;\rightarrow\; \overleftarrow{\mathbf{P}_{r(s)} \;\rightarrow\; \mathbf{S}_{rs}} \;\rightarrow\; \mathbf{P}^p_{r(s)} \;\rightarrow\; \mathbf{S}^p_{(r)}\f$ |

    Rows 1-6 describe simulation modes that include primary emission only, without or with media.
    The last mode in this list (row 6) includes iteration over the DMS during primary emission.
    Rows 7-11 describe modes that include dust and/or gas emission, possibly iterating over the RF
    during secondary emission to handle DSE (rows 9-11). If one or more of the dust media have a
    DMS, there are modes supporting iteration over primary emission (rows 8 and 10) or even over
    both primary and secondary emission (row 11).

    The simulation modes in rows 7-11 support medium components with SDMS by updating the medium
    state for these components at the end of the relevant segments. It is meaningful to update an
    SDMS at the end of a segment that performs peel-off to the instruments because the update does
    not affect the outcome of a primary emission segment.

    */
class MonteCarloSimulation : public Simulation
{
    /** The enumeration type indicating the simulation mode, which determines the overall structure
        of the simulation and its capabilities. The choice made for the simulation mode has a
        significant impact on which options are allowed or required in the simulation's
        configuration. See the class header for more information. */
    ENUM_DEF(SimulationMode, OligoNoMedium, OligoExtinctionOnly, NoMedium, ExtinctionOnly, LyaExtinctionOnly,
             DustEmission, GasEmission, DustAndGasEmission)
        ENUM_VAL(SimulationMode, OligoNoMedium, "No medium - oligochromatic regime (a few discrete wavelengths)")
        ENUM_VAL(SimulationMode, OligoExtinctionOnly,
                 "Extinction only - oligochromatic regime (a few discrete wavelengths)")
        ENUM_VAL(SimulationMode, NoMedium, "No medium (primary sources only)")
        ENUM_VAL(SimulationMode, ExtinctionOnly, "Extinction only (no secondary emission)")
        ENUM_VAL(SimulationMode, LyaExtinctionOnly, "Extinction only with Lyman-alpha line transfer")
        ENUM_VAL(SimulationMode, DustEmission, "With secondary emission from dust")
        ENUM_VAL(SimulationMode, GasEmission, "With secondary emission from gas")
        ENUM_VAL(SimulationMode, DustAndGasEmission, "With secondary emission from dust and gas")
    ENUM_END()

    ITEM_CONCRETE(MonteCarloSimulation, Simulation, "a Monte Carlo simulation")

        PROPERTY_ENUM(simulationMode, SimulationMode, "the overall simulation mode")
        ATTRIBUTE_DEFAULT_VALUE(simulationMode, "ExtinctionOnly")
        ATTRIBUTE_INSERT(
            simulationMode,
            "simulationModeOligoNoMedium:Oligochromatic,NoMedium;"
            "simulationModeOligoExtinctionOnly:Oligochromatic,ExtinctionOnly;"
            "simulationModeNoMedium:Panchromatic,NoMedium;"
            "simulationModeExtinctionOnly:Panchromatic,ExtinctionOnly;"
            "simulationModeLyaExtinctionOnly:Lya,Panchromatic,ExtinctionOnly;"
            "simulationModeDustEmission:Panchromatic,DustEmission,Emission,RadiationField;"
            "simulationModeGasEmission:Panchromatic,GasEmission,Emission,RadiationField;"
            "simulationModeDustAndGasEmission:Panchromatic,DustEmission,GasEmission,Emission,RadiationField")

        PROPERTY_BOOL(iterateMediumState, "iterate to self-consistently determine the medium state")
        ATTRIBUTE_DEFAULT_VALUE(iterateMediumState, "false")
        ATTRIBUTE_RELEVANT_IF(iterateMediumState, "simulationModeExtinctionOnly|Emission")
        ATTRIBUTE_DISPLAYED_IF(iterateMediumState, "Level3")
        ATTRIBUTE_INSERT(iterateMediumState,
                         "(simulationModeExtinctionOnly|Emission)&iterateMediumState:DynamicState,RadiationField")

        PROPERTY_BOOL(iterateSecondaryEmission, "iterate to self-consistently determine secondary emission")
        ATTRIBUTE_DEFAULT_VALUE(iterateSecondaryEmission, "false")
        ATTRIBUTE_RELEVANT_IF(iterateSecondaryEmission, "Emission")
        ATTRIBUTE_DISPLAYED_IF(iterateSecondaryEmission, "Level2")
        ATTRIBUTE_INSERT(iterateSecondaryEmission, "Emission&iterateSecondaryEmission:DynamicEmission")

        PROPERTY_DOUBLE(numPackets, "the default number of photon packets launched per simulation segment")
        ATTRIBUTE_MIN_VALUE(numPackets, "[0")
        ATTRIBUTE_MAX_VALUE(numPackets, "1e19]")
        ATTRIBUTE_DEFAULT_VALUE(numPackets, "1e6")

        PROPERTY_ITEM(cosmology, Cosmology, "the cosmology parameters")
        ATTRIBUTE_DEFAULT_VALUE(cosmology, "LocalUniverseCosmology")
        ATTRIBUTE_DISPLAYED_IF(cosmology, "Level2")

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

    /** This function constructs a SecondarySourceSystem object if the simulation configuration
        requires secondary emission. */
    void setupSelfAfter() override;

    //======== Getters for Non-Discoverable Properties =======

public:
    /** Returns the Configuration object for this simulation hierarchy. */
    Configuration* config() const;

    //======================== Other Functions =======================

protected:
    /** This function actually runs the simulation, assuming setup has been completed. It
        determines the appropriate execution flow from the simulation's configuration (see the
        table in the class header documentation) and invokes other functions in this class to
        consecutively perform all segments of the simulation, including primary source segments,
        secondary source segments, and any required iterations. For each segment, these invoked
        functions will tell the source system to prepare for launching photon packets (in serial
        code), and then cause photon packets to be launched and traced through their life cycles
        (in appropriately parallelized code depending on the run-time environment and the
        command-line options). */
    void runSimulation() override;

private:
    /** This function runs a final primary source emission segment including peel-off towards the
        instruments. It records radiation field contributions if the configuration requires it
        (e.g., because secondary emission must be calculated, or because the user configured probes
        to directly output radiation field information). If the configuration also includes
        emission, any semi-dynamic state media in the simulation receive an update request at the
        end of the segment.

        Using the notation of the table in the class header documentation, this function implements
        the execution flow \f$\mathbf{P}^p\f$, \f$\mathbf{P}^p_{(r)}\f$, or
        \f$\mathbf{P}^p_{r(s)}\f$, depending on the relevant configuration options. */
    void runPrimaryEmission();

    /** This function runs a final secondary source emission segment including peel-off towards the
        instruments. It records radiation field contributions if the configuration requires it
        (e.g., because the user configured probes to directly output radiation field information).

        Using the notation of the table in the class header documentation, this function implements
        the \f$\mathbf{S}^p_{(r)}\f$ execution flow. */
    void runSecondaryEmission();

    /** This function iteratively runs consecutive primary source emission segments to
        self-consistently calculate the radiation field taking into account one or more recipes for
        updating the dynamic medium state as a function of the radiation field. See the
        DynamicStateRecipe class for more information.

        Each segment in the iteration tracks photon packets through the medium and subsequently
        updates the medium state based on the newly established radiation field. The number of
        photon packets launched for each segment can be configured as a multiplication factor on
        the default number of primary photon packets. The minimum and maximum number of iterations
        can also be specified as configuration options. Within these limits, the actual number of
        iterations performed is determined by convergence criteria defined by the configured
        dynamic medium state recipe(s).

        The photon packets emitted by this function do not perform peel-off to the instruments
        because the radiation field has not yet converged. After this function returns, a final
        segment of primary emission with peel-off must be performed by calling the
        runPrimaryEmission() function.

        Using the notation of the table in the class header documentation, this function implements
        the \f$\overleftarrow{\mathbf{P}_{rs}}\f$ execution flow. */
    void runPrimaryEmissionIterations();

    /** This function iteratively runs consecutive secondary source emission segments to
        self-consistently calculate the radiation field in models where the secondary emission
        spectrum and the total (primary plus secondary) radiation field significantly depend on
        each other. A typical example is a model where the dust medium is sufficiently opaque at
        infrared wavelengths to cause significant self-absorption and thus additional heating.

        Each segment in the iteration tracks photon packets through the medium to calculate the
        radiation field resulting from the secondary sources in the simulation. The subsequent
        segment then recalculates the secondary emission spectrum based on the newly established
        total (primary and secondary) radiation field. The number of photon packets launched for
        each segment can be configured as a multiplication factor on the default number of primary
        photon packets. The minimum and maximum number of iterations can also be specified as
        configuration options. Within these limits, the actual number of iterations performed is
        determined by user-configured convergence criteria. Specifically, convergence is reached
        when (a) the absorbed secondary luminosity is less than a given fraction of the absorbed
        primary luminosity, \em OR (b) the absorbed secondary luminosity has changed by less than a
        given fraction compared to the previous iteration.

        \note In the current implementation, the convergence criteria are based solely on dust
        absorption. For proper operation, this function thus requires the simulated model to
        contain one or more dust components and have dust emission enabled.

        Any semi-dynamic state media in the simulation receive an update request at the end of each
        segment.

        The photon packets emitted by this function do not perform peel-off to the instruments
        because the radiation field has not yet converged. After this function returns, a final
        segment of secondary emission with peel-off must be performed by calling the
        runSecondaryEmission() function.

        Using the notation of the table in the class header documentation, this function implements
        the \f$\overleftarrow{\mathbf{S}_{r(s)}}\f$ execution flow. */
    void runSecondaryEmissionIterations();

    /** This function iteratively runs consecutive primary and secondary source emission segments
        to self-consistently calculate the radiation field in models that have both a dynamic
        medium state (the extinction properties depend on the radiation field) and dynamic
        secondary emission (the emission spectrum depends on the radiation field), and the user
        requested to include the primary emission in secondary emission iterations. This would be
        meaningful in case the changes in the radiation field caused by secondary emission could
        sufficiently influence the medium state to significantly affect the results of primary
        emission.

        The iteration is considered to be converged only when both the criteria for primary
        emission iteration and those for secondary emission iteration are satisfied; for more
        information see the runPrimaryEmissionIterations() and runSecondaryEmissionIterations()
        functions.

        Any semi-dynamic state media in the simulation receive an update request at the end of each
        primary segment, and the dynamic state recipes are requested to update the dynamic state
        after each secondary emission segment.

        The photon packets emitted by this function do not perform peel-off to the instruments
        because the radiation field has not yet converged. After this function returns, a final
        segment of secondary emission with peel-off must be performed by calling the
        runSecondaryEmission() function.

        Using the notation of the table in the class header documentation, this function implements
        the \f$\overleftarrow{\mathbf{P}_{r(s)} \;\rightarrow\; \mathbf{S}_{rs}}\f$ execution flow.
        */
    void runMergedEmissionIterations();

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
        probability that it is emitted in any other direction. If the source has a nonzero
        velocity, the wavelength of the peel-off photon packet is Doppler-shifted for the new
        direction.

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
        perceived wavelength in the cell under consideration, taking into account the velocity of
        the medium in that cell.

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
