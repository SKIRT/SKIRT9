/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NONLTELINEGASMIX_H
#define NONLTELINEGASMIX_H

#include "EmittingGasMix.hpp"

////////////////////////////////////////////////////////////////////

/** The NonLTELineGasMix class describes the material properties related to selected transitions in
    selected molecules and atoms. For each supported species, the current implementation includes a
    number of rotational energy levels (quantum number \f$J\f$) at the base vibrational level
    (quantum number \f$v=0\f$) for molecules and electronic energy levels and hyperfine-split
    levels for atoms, and supports the allowed transitions between these levels. Vibrational energy
    levels and the corresponding rovibrational transitions may be added later. The class properties
    allow configuring the species and the number of transitions to be considered.

    For each supported transition, the emission luminosity and absorption opacity in a given cell
    are determined from the gas properties defined in the input model and the local radiation field
    calculated by the simulation. The class performs an iterative, non-LTE calculation including
    the effects of collisional transitions (excitation and de-excitation with one or more types of
    interaction partners) and photonic transitions (spontaneous emission, absorption, and induced
    emission). This allows establishing the energy level distribution (population) for a wide range
    of material densities and radiation fields.


    <b>Supported species and transitions</b>

    The current implementation supports the following molecular or atomic species:

    <b>Molecules</b>
    The properties relevant to the calculation of the level populations and collision partners 
    of the following molecules are based on the LAMDA database (
    https://home.strw.leidenuniv.nl/~moldata/; Schöier et al. 2005; van der Tak et al. 2020).
    According to the user's choice, the number of energy levels can be set by the 
    \em maxEnergyLevels property within the range of supported energy levels for each species.

    - \c Two-level \c test molecule (TT): a fictive test molecule (called TT for our purposes)
    that has just two rotational energy levels with a corresponding transition line at 1666.67
    \f$\mu\mathrm{m}\f$. The single collisional interaction partner is molecular hydrogen. The
    properties of this molecule are defined by van Zadelhoff et al. 2002 for use with the first
    benchmark problem described there.

    - \c Hydroxyl radical (OH): includes rotational energy levels up to \f$J=25/2\f$ without
    hyperfine-split levels. The corresponding transition lines are at wavelengths from 10
    \f$\mu\mathrm{m}\f$ to 2 \f$\mathrm{m}\f$. The collisional interaction partners include
    molecular and neutral atomic hydrogen and neutral atomic Helium.

    - \c Hydroxyl radical (OHhfs): includes rotational energy levels up to \f$J=5/2\f$ including
    hyperfine-split levels. The corresponding transition lines are at wavelengths from 34.6
    \f$\mu\mathrm{m}\f$ to 186 \f$\mathrm{mm}\f$. The single collisional interaction partner is
    molecular hydrogen.

    - \c Formyl cation (HCO+): includes rotational energy levels up to \f$J=21\f$. The
    corresponding transition lines are at wavelengths from 160.3 to 3361 \f$\mu\mathrm{m}\f$. The
    single collisional interaction partner is molecular hydrogen.

    - \c HydrogenCyanide cation (HCN): includes rotational energy levels up to \f$J=25\f$. The
    corresponding transition lines are at wavelengths from 135.6 to 3382 \f$\mu\mathrm{m}\f$. The
    collisional interaction partners include molecular hydrogen and electrons.

    - \c Carbon \c monoxide (CO): includes rotational energy levels up to \f$J=40\f$. The
    corresponding transition lines are at wavelengths from 65.7 to 2601 \f$\mu\mathrm{m}\f$. The
    single collisional interaction partner is molecular hydrogen.

    - \c Molecular \c hydrogen (H2): includes rotational energy levels up to \f$J=31\f$. The
    corresponding transition lines are at wavelengths from 28.2 to 3.3 \f$\mu\mathrm{m}\f$. The
    collisional interaction partners include molecular hydrogen, neutral atomic hydrogen, ionized
    atomic hydrogen, and helium. 
    
    <b> An assumption of molecular hydrogen:</b> Molecular hydrogen exists in ortho (J is odd) and
    para (J is even) flavors, and these flavors have different collisional coefficients.
    The current implementation assumes that the ratio of ortho to para is 3 to 1.

    <b>Atoms</b>
    Atomic line emission from C, O, N, S, Si, Ne, Na, Mg, Fe, and Ar in ionization states ranging from
    neutral up to ten times ionized is included. The atomic data required for calculating level 
    populations and electron collisional transition rates are adopted from the CHIANTI database
    (https://www.chiantidatabase.org/chianti_direct_access.html; Dere et al. 1997; Del Zanna et al. 2021).

    In addition, collisional transition rates with other collision partners are available for CI, CII, 
    OI, OIII, NII, SiI, and SI, based on the LAMDA database. The available collision partners are neutral, 
    ionized, and molecular hydrogen, as well as neutral helium for CI and OI; atomic and molecular 
    hydrogen for CII; and atomic hydrogen for OIII, NII, SiI, and SI.

    Level populations can be calculated from the ground state up to the 20th excited level, including all
    hyperfine-split sublevels within each electronic state. The corresponding radiative transitions span
    wavelengths from the ultraviolet to the submillimeter regime. The maximum number of energy levels
    included in the calculation can be specified by the user through the \em maxEnergyLevels property,
    up to this limit. For example, users may be interested in only a famous following far-infrared line
    emission rather than UV and optical lines: CII (158 \f$\mu\mathrm{m}\f$), 
    OI (63 and 145 \f$\mu\mathrm{m}\f$), NII (122 and 205 \f$\mu\mathrm{m}\f$),
    CI (370 and 609 \f$\mu\mathrm{m}\f$), OIII (52 and 88 \f$\mu\mathrm{m}\f$). 
    \em maxEnergyLevels should be set to 2 for CII, 3 for OI, NII, CI, and OIII.

    <b> Note on highly ionized species: </b> The current implementation does not include photoionization
    and recombination processes. Consequently, the level populations of highly ionized species may not be
    accurately determined in environments with strong ionizing radiation fields.

    <b> Neglecting irrelevant transitions:</b> Transitions that have negligible impact on the level
    populations can be safely ignored to improve computational efficiency by setting the 
    \em lowestBranchingRatio property to a value greater than zero. The default value of 0.0 
    includes all transitions, while a value of 0.001, for example, would exclude the calculation
    of emissivity and opacity of the transitions that have a branching ratio below 0.1% of the
    total decay rate from the upper level. This can reduce the number of transitions considered in
    the calculations and avoide unnecessary computations for transitions that have a negligible
    effect on the level populations.

    <b>Configuring the simulation</b>

    Simulations that include gas represented by the NonLTELineGasMix often also include dust,
    although this is not a requirement. In any case, the simulation must have one or more primary
    sources that trigger the molecular lines directly (e.g. the cosmic microwave background) or
    indirectly (e.g. by heating the dust and thus causing thermal dust emission), or both. The \em
    simulationMode should be set to "GasEmission" or "DustAndGasEmission" correspondingly.

    Because the secondary emission and the radiation field are calculated self-consistently, \em
    iterateSecondaryEmission should always be enabled. If the wavelength range of the primary
    source(s) significantly overlaps the emission lines of the species under consideration (in
    other words, if the opacity at the line wavelengths significantly affects the primary radiation
    field), then \em includePrimaryEmission in IterationOptions must be enabled as well. If not,
    this property can be left disabled. In both cases, \em iteratePrimaryEmission can be left
    disabled unless iteration is required for other media (for example, to self-consistently
    calculate radiative dust destruction).

    The radiation field wavelength grid should properly resolve the UV, optical, and infrared
    wavelength range (in case the simulation includes dust) and the wavelength ranges of the
    supported transition lines. Separate instruments can be configured for the relevant wavelength
    ranges, e.g. using a logarithmic grid for the continuum spectrum and linear grids for the line
    profiles.

    During the initial primary emission segment, the simulation determines the radiation field
    resulting from the primary sources, dust attenuation (if present) and molecular line absorption
    based on default equilibrium level populations. The resulting radiation field allows a first
    estimation of the level populations for each spatial grid cell; this calculation happens in the
    updateSpecificState() function.

    During secondary emission, the simulation takes into account emission from all configured media
    in addition to absorption by these same media, including molecular lines in both cases. The
    previously stored level populations allow calculating the line emission spectrum and the line
    absorption cross sections. This results in an updated radiation field, which will in turn
    influence the level populations (and for high optical depths, possibly the dust temperature),
    which in turn influences the secondary emission spectra. In order to obtain a self-consistent
    result, the simulation must therefore iterate over secondary emission.

    The input model must provide values for the spatial distribution of several medium properties,
    including the number density of the species under consideration, the number density of any
    relevant collisional partner species, the kinetic gas temperature, and the turbulence velocity.
    These values remain constant during the simulation. Most often, this information will be read
    from an input file by associating the NonLTELineGasMix with a subclass of ImportedMedium.
    For that medium component, the ski file attribute \em importTemperature <b>must</b> be set to
    'true', and \em importMetallicity and \em importVariableMixParams must be left at 'false'. The
    additional columns required by the material mix are automatically imported and are expected
    <b>after</b> all other columns. For example, if bulk velocities are also imported for this
    medium component (i.e. \em importVelocity is 'true'), the column order would be \f[ ...,
    n_\mathrm{mol}, T_\mathrm{kin}, v_\mathrm{x}, v_\mathrm{y}, v_\mathrm{z}, n_\mathrm{col1} [,
    n_\mathrm{col1}, ...], v_\mathrm{turb}\f]

    For basic testing purposes, the NonLTELineGasMix can also be associated with a geometric
    medium component. The geometry then defines the spatial density distribution of the species
    under consideration (i.e. \f$n_\mathrm{mol}\f$), and the NonLTELineGasMix configuration
    properties specify a fixed default value for the other properties that will be used across the
    spatial domain. In this case, the number densities of the collisional partners are defined by a
    constant multiplier relative to \f$n_\mathrm{mol}\f$.

    Although the NonLTELineGasMix couples the level populations to the radiation field, in some cases
    the radiation field does not significantly affect the level populations (e.g., UV and optical 
    lines from atomic species). In such cases, the user can disable updating the level populations by
    setting the \em updateDynamicStatesFlag property to 'false'.

    <b>Thermal motion and turbulence</b>

    The thermal velocity in a medium of particles with mass \f$m\f$ at temperature
    \f$T_\mathrm{kin}\f$ is defined as \f[ v_\mathrm{th} = \sqrt{\frac{2kT_\mathrm{kin}}{m}}. \f]
    This value corresponds to the most probable particle speed, i.e. the point where the
    probability distribution of the velocity vector norm reaches its maximum value. One often
    considers an additional source of line broadening caused by subgrid processes other than those
    corresponding to the kinetic temperature. This motion is characterized by the turbulent
    velocity \f$v_\mathrm{turb}\f$. Assuming a Gaussian line profile, the total velocity dispersion
    \f$v_\mathrm{s}\f$ (the standard deviation of the Gaussian in velocity space) is then defined
    through \f[ \sqrt{2}\,v_\mathrm{s} = \sqrt{v_\mathrm{th}^2 + v_\mathrm{turb}^2}. \f] We can
    artificially combine the effect of both thermal motion and turbulence into an effective
    temperature, \f[ T_\mathrm{eff} = T_\mathrm{kin} + \frac{m v_\mathrm{turb}^2}{2k}. \f]
    If the Gaussian profile is not expressed when convolving the radiation field with the Gaussian
    profile, the simulation may fail and output an error. Sometimes, the emission is not critical 
    for estimating the level populations. In this case, the user can choose to take a warning instead
    of an error by setting the \em errorForGaussianIntegral property to 'false'.

    <b>Numerical convergence</b>

    As mentioned above, a simulation including this material mix performs iterations to
    self-consistently calculate the radiation field and the medium state (i.e. the energy level
    populations that drive the line emission). This class offers three user-configurable properties
    that specify the convergence criteria for this iterative process.

    The first two properties configure a criterion based on statistics per spatial cell. The \em
    maxChangeInLevelPopulations property specifies the maximum relative change between consecutive
    iterations in the level populations for a given spatial cell for that cell to be considered
    converged. The \em maxFractionNotConvergedCells property then specifies the maximum fraction of
    cells (relative to the total number of cells in the simulation) that may be left not converged
    for the whole spatial domain to be considered converged.

    The third property configures a global criterion. The \em maxChangeInGlobalLevelPopulations
    property specifies the maximum relative change between consecutive iterations in the global
    level populations, accumulated over the complete spatial domain, for the spatial domain to be
    considered converged. In other words, convergence has been reached if \f[ \epsilon =
    \underset{i}{\mathrm{max}} \left| \frac{N_i^{(k)} - N_i^{(k-1)}} {N_i^{(k-1)}} \right| \le
    \epsilon_\mathrm{max}\f] where \f$N_i^{(k)}\f$ denotes the total population across the spatial
    domain for energy level \f$i\f$ at iteration \f$k\f$ and \f$\epsilon_\mathrm{max}\f$ is the
    user-configured value of the \em maxChangeInGlobalLevelPopulations property.

    Both criteria must be satisfied simultaneously. Howver, the user can effectively disable a
    criterion by specifying very liberal maximum values.

    <b>Level populations</b>

    We denote the populations for the \f$N\f$ supported energy levels as \f$n_i\f$, with indices
    \f$i=0,N-1\f$ increasing from lowest to highest energy state, and with \f$\sum_i n_i =
    n_\mathrm{mol}\f$. The level populations form the basis to calculate the emission luminosities
    and the absorption opacities at later stage. Their values are obtained by solving the set of
    statistical equilibrium equations given by

    \f[ \sum_{j>i} \Big[n_jA_{ji} + (n_jB_{ji} - n_iB_{ij})J_{\lambda,ji})\Big]
    - \sum_{j<i} \Big[n_iA_{ij} + (n_iB_{ij} - n_jB_{ji})J_{\lambda,ij})\Big]
    + \sum_{j\neq i} \Big[n_jC_{ji} - n_iC_{ij}\Big]=0, \quad i=0,N-1\f]

    where \f$J_{\lambda,ul}\f$ is the mean radiation field intensity integrated over the line
    profile corresponding to the transition from upper energy level \f$u\f$ to lower energy level
    \f$l\f$. Furthermore, \f$A_{ul}\f$, \f$B_{ul}\f$, \f$B_{lu}\f$ are the Einstein coefficients
    for spontaneous emission, induced emission, and absorption, and \f$C_{ul}\f$, \f$C_{lu}\f$ the
    coefficients for collisional excitation and de-excitation. These coefficients are taken from
    the literature.

    The Einstein coefficients \f$A_{ul}\f$, \f$B_{ul}\f$, and \f$B_{lu}\f$ are constant for each
    transition, and the \f$B_{ul}\f$ and \f$B_{lu}\f$ coefficients can be obtained from the
    \f$A_{ul}\f$ coefficients through \f[ B_{ul} = \frac{\lambda_{ul}^5}{2 hc^2}A_{ul} \f] and \f[
    B_{lu} = \frac{g_u}{g_l} B_{ul} = \frac{g_u}{g_l} \frac{\lambda_{ul}^5}{2 hc^2}A_{ul}, \f]
    where \f$g_u\f$ and \f$g_l\f$ represent the degeneracy of the upper and lower energy levels,
    respectively. The collisional coefficients \f$C_{ul}\f$ and \f$C_{lu}\f$ depend on the number
    density of the collisional partner and on the kinetic temperature \f$T_\mathrm{kin}\f$ of the
    gas. These coefficients are related by \f[ \frac{C_{lu}}{C_{ul}} = \frac{g_u}{g_l}
    \exp(-E_{ul}/kT_\mathrm{kin}).\f]

    <b>Emission</b>

    The integrated line luminosity \f$L_{ul}\f$ corresponding to the transition from upper energy
    level \f$u\f$ to lower energy level \f$l\f$ for a given spatial cell is given by

    \f[ L_{ul} = \frac{hc}{\lambda_{ul}} A_{ul} n_u V_\mathrm{cell}, \f]

    where \f$\lambda_{ul}\f$ is the transition wavelength and \f$V_\mathrm{cell}\f$ is the volume
    of the cell. The SKIRT framework automatically adds a Gaussian line profile assuming a thermal
    velocity corresponding to the effective temperature in the cell and the mass of a molecule, in
    addition to any Doppler shifts caused by the bulk velocity in the cell.

    <b>Absorption</b>

    The absorption opacity \f$k^\text{abs}(\lambda)\f$ as a function of wavelength \f$\lambda\f$ is
    given by

    \f[ k^\text{abs}(\lambda) = \sum_{ul} \, \frac{h c}{4\pi \lambda_{ul}}
    (n_l\,B_{lu}-n_u\,B_{ul}) \,\phi_{ul}(\lambda) \f]

    where the sum runs over all supported transitions, \f$u\f$ and \f$l\f$ are the indices of the
    energy levels before and after the transition, and \f$\phi_{ul}(\lambda)\f$ is the Gaussian
    profile of the line corresponding to the transition.

    The calculation explicitly includes the Gaussian line profile caused by thermal motion of the
    molecules in the local bulk velocity frame of the cell. Moreover, the absorption profiles of
    all supported lines are superposed on top of each other. In practice, the implementation
    includes just the terms that have a significant contribution at any given wavelength.

    <b>Limiting negative optical depth</b>

    The formula for the absorption opacity given in the previous section yields a negative value in
    case the medium exhibits stimulated emission. The explicit absorption technique (see the
    MonteCarloSimulation class) supports the corresponding negative optical depths along a photon
    path. However, numerical problems arise if the negative optical depth magnitude becomes too
    high. This may happen, for example, as a result of Monte Carlo noise, or in the early stages
    when a simulation has not yet converged.

    Therefore, this class imposes a user-configurable lower limit on the (negative) optical depth
    along a cell diagonal before returning an absorption opacity. The default limit is -2.

    <b>Storing mean intensities</b>

    As discussed above, the mean radiation field intensity \f$J_{\lambda,ul}\f$ is determined for
    each transition line as a prerequisite for calculating the level populations in a given cell.
    These values are obtained by integrating the radiation field, as discretized on the
    simulation's radiation field wavelength grid, over the appropriate Gaussian line profile. There
    is no need to keep these values around for the regular operation of the code. However, they can
    form a relevant diagnostic to assess the level of Monte Carlo noise in the simulation or to
    evaluate iterative convergence. Therefore, if the \em storeMeanIntensities flag is turned on,
    the mean radiation field intensity at each transition line is stored in the medium state for
    each cell, so that the information can be probed using the CustomStateProbe class.

    <b>Providing initial level populations</b>
    The \em initialLevelPopsCase property specifies how the initial level populations are assigned.
    By default, the populations in each cell are initialized assuming a Boltzmann distribution at
    the local gas temperature (i.e., LTE conditions; \em initialLevelPopsCase="LTE"). This
    assumption may significantly differ from the final non-LTE solution, although it often provides
    a convenient starting point that accelerates convergence. For example, many optical and UV 
    emission lines from the ionized gas originate from excited levels whose populations are 
    determined by the balance between collisional transitions and spontaneous radiative decay. In 
    such cases, the user can initialize the populations using the collisionally excited solution
    by specifying \em initialLevelPopsCase="CollisionallyExcited". Alternatively, the users can 
    set custom initial level populations using a text file by specifying 
    \em initialLevelPopsCase="Custom". This option is particularly useful when restarting a SKIRT
    simulation from a previously converged solution.  If the \em initialLevelPopsFilename string is
    nonempty, it specifies the name of a text column file with a column for each energy level and
    a row for each spatial cell in the simulation. Specifically, the first column lists a cell 
    index that is ignored, and remaining columns list the relative population for each energy level
    in units of number density (the values are scaled to the total number density in the cell, so
    the specific units don't really matter). The rows must exactly match the number and ordering
    of the cells in the simulation's spatial grid (see below).

    This input format is designed such that the text file can easily be produced by a
    CustomStateProbe instance in a previous simulation. For example, assuming 9 energy levels, one
    could configure a probe as follows:

        <CustomStateProbe probeName="populations" indices="2-10" probeAfter="Run">
            <form type="Form">
                <PerCellForm/>
            </form>
        </CustomStateProbe>

    The requirement to have an identical spatial grid in both simulations is easily met for 1D, 2D
    and 3D grids with a regular mesh in each spatial direction. For octree and binary tree grids
    the precise cell hierarchy is usually determined based on random samples of the medium density
    distribution, resulting in a slightly different structure for each run. To work around this
    problem, one needs to output the topology of the hierarchical grid in the first simulation
    using a TreeSpatialGridTopologyProbe instance, and load this topology in subsequent simulations
    using a FileTreeSpatialGrid instance.

    */
class NonLTELineGasMix : public EmittingGasMix
{
    /** The enumeration type indicating the molecular or atomic species represented by a given
        NonLTELineGasMix instance. See the class header for more information. */
    ENUM_DEF(Species, Test, Hydroxyl, HydroxylHFS, Formyl, HydrogenCyanide, CarbonMonoxide, MolecularHydrogen,
             AtomicCarbon, IonizedCarbon, DoublyIonizedCarbon, TriplyIonizedCarbon, FourTimesIonizedCarbon,
             FiveTimesIonizedCarbon, AtomicNitrogen, IonizedNitrogen, DoublyIonizedNitrogen, TriplyIonizedNitrogen,
             FourTimesIonizedNitrogen, FiveTimesIonizedNitrogen, SixTimesIonizedNitrogen, AtomicOxygen, IonizedOxygen,
             DoublyIonizedOxygen, TriplyIonizedOxygen, FourTimesIonizedOxygen, FiveTimesIonizedOxygen,
             SixTimesIonizedOxygen, SevenTimesIonizedOxygen, AtomicNeon, IonizedNeon, DoublyIonizedNeon,
             TriplyIonizedNeon, FourTimesIonizedNeon, FiveTimesIonizedNeon, SixTimesIonizedNeon, SevenTimesIonizedNeon,
             EightTimesIonizedNeon, NineTimesIonizedNeon, IonizedSodium, DoublyIonizedSodium, TriplyIonizedSodium,
             FourTimesIonizedSodium, FiveTimesIonizedSodium, SixTimesIonizedSodium, SevenTimesIonizedSodium,
             EightTimesIonizedSodium, NineTimesIonizedSodium, TenTimesIonizedSodium, IonizedMagnesium,
             DoublyIonizedMagnesium, TriplyIonizedMagnesium, FourTimesIonizedMagnesium, FiveTimesIonizedMagnesium,
             SixTimesIonizedMagnesium, SevenTimesIonizedMagnesium, EightTimesIonizedMagnesium,
             NineTimesIonizedMagnesium, TenTimesIonizedMagnesium, AtomicSilicon, IonizedSilicon, DoublyIonizedSilicon,
             TriplyIonizedSilicon, FourTimesIonizedSilicon, FiveTimesIonizedSilicon, SixTimesIonizedSilicon,
             SevenTimesIonizedSilicon, EightTimesIonizedSilicon, NineTimesIonizedSilicon, TenTimesIonizedSilicon,
             AtomicSulfur, IonizedSulfur, DoublyIonizedSulfur, TriplyIonizedSulfur, FourTimesIonizedSulfur,
             FiveTimesIonizedSulfur, SixTimesIonizedSulfur, SevenTimesIonizedSulfur, EightTimesIonizedSulfur,
             NineTimesIonizedSulfur, TenTimesIonizedSulfur, DoublyIonizedArgon, TriplyIonizedArgon,
             FourTimesIonizedArgon, FiveTimesIonizedArgon, SixTimesIonizedArgon, SevenTimesIonizedArgon,
             EightTimesIonizedArgon, NineTimesIonizedArgon, TenTimesIonizedArgon, IonizedIron, DoublyIonizedIron,
             TriplyIonizedIron, FourTimesIonizedIron, FiveTimesIonizedIron, SixTimesIonizedIron, SevenTimesIonizedIron,
             EightTimesIonizedIron, NineTimesIonizedIron, TenTimesIonizedIron)
        ENUM_VAL(Species, Test, "Fictive two-level test molecule (TT)")
        ENUM_VAL(Species, Hydroxyl, "Hydroxyl radical (OH)")
        ENUM_VAL(Species, HydroxylHFS, "Hydroxyl radical (OH) with hyperfine structure")
        ENUM_VAL(Species, Formyl, "Formyl cation (HCO+)")
        ENUM_VAL(Species, HydrogenCyanide, "Hydrogen cyanide (HCN)")
        ENUM_VAL(Species, CarbonMonoxide, "Carbon monoxide (CO)")
        ENUM_VAL(Species, MolecularHydrogen, "Molecular hydrogen (H2)")

        ENUM_VAL(Species, AtomicCarbon, "Atomic carbon (CI)")
        ENUM_VAL(Species, IonizedCarbon, "Ionized carbon (CII)")
        ENUM_VAL(Species, DoublyIonizedCarbon, "Doubly Ionized Carbon (CIII)")
        ENUM_VAL(Species, TriplyIonizedCarbon, "Triply Ionized Carbon (CIV)")
        ENUM_VAL(Species, FourTimesIonizedCarbon, "Four times Ionized Carbon (CV)")
        ENUM_VAL(Species, FiveTimesIonizedCarbon, "Five times Ionized Carbon (CVI)")

        ENUM_VAL(Species, AtomicNitrogen, "Atomic nitrogen (NI)")
        ENUM_VAL(Species, IonizedNitrogen, "Ionized nitrogen (NII)")
        ENUM_VAL(Species, DoublyIonizedNitrogen, "Doubly Ionized Nitrogen (NIII)")
        ENUM_VAL(Species, TriplyIonizedNitrogen, "Triply Ionized Nitrogen (NIV)")
        ENUM_VAL(Species, FourTimesIonizedNitrogen, "Four times Ionized Nitrogen (NV)")
        ENUM_VAL(Species, FiveTimesIonizedNitrogen, "Five times Ionized Nitrogen (NVI)")
        ENUM_VAL(Species, SixTimesIonizedNitrogen, "Six times Ionized Nitrogen (NVII)")

        ENUM_VAL(Species, AtomicOxygen, "Atomic oxygen (OI)")
        ENUM_VAL(Species, IonizedOxygen, "Ionized oxygen (OII)")
        ENUM_VAL(Species, DoublyIonizedOxygen, "Doubly Ionized Oxygen (OIII)")
        ENUM_VAL(Species, TriplyIonizedOxygen, "Triply Ionized Oxygen (OIV)")
        ENUM_VAL(Species, FourTimesIonizedOxygen, "Four times Ionized Oxygen (OV)")
        ENUM_VAL(Species, FiveTimesIonizedOxygen, "Five times Ionized Oxygen (OVI)")
        ENUM_VAL(Species, SixTimesIonizedOxygen, "Six times Ionized Oxygen (OVII)")
        ENUM_VAL(Species, SevenTimesIonizedOxygen, "Seven times Ionized Oxygen (OVIII)")

        ENUM_VAL(Species, AtomicNeon, "Atomic Neon (NeI)")
        ENUM_VAL(Species, IonizedNeon, "Ionized Neon (NeII)")
        ENUM_VAL(Species, DoublyIonizedNeon, "Doubly Ionized Neon (NeIII)")
        ENUM_VAL(Species, TriplyIonizedNeon, "Triply Ionized Neon (NeIV)")
        ENUM_VAL(Species, FourTimesIonizedNeon, "Four times Ionized Neon (NeV)")
        ENUM_VAL(Species, FiveTimesIonizedNeon, "Five times Ionized Neon (NeVI)")
        ENUM_VAL(Species, SixTimesIonizedNeon, "Six times Ionized Neon (NeVII)")
        ENUM_VAL(Species, SevenTimesIonizedNeon, "Seven times Ionized Neon (NeVIII)")
        ENUM_VAL(Species, EightTimesIonizedNeon, "Eight times Ionized Neon (NeIX)")
        ENUM_VAL(Species, NineTimesIonizedNeon, "Nine times Ionized Neon (NeX)")

        ENUM_VAL(Species, IonizedSodium, "Ionized Sodium (NaII)")
        ENUM_VAL(Species, DoublyIonizedSodium, "Doubly Ionized Sodium (NaIII)")
        ENUM_VAL(Species, TriplyIonizedSodium, "Triply Ionized Sodium (NaIV)")
        ENUM_VAL(Species, FourTimesIonizedSodium, "Four times Ionized Sodium (NaV)")
        ENUM_VAL(Species, FiveTimesIonizedSodium, "Five times Ionized Sodium (NaVI)")
        ENUM_VAL(Species, SixTimesIonizedSodium, "Six times Ionized Sodium (NaVII)")
        ENUM_VAL(Species, SevenTimesIonizedSodium, "Seven times Ionized Sodium (NaVIII)")
        ENUM_VAL(Species, EightTimesIonizedSodium, "Eight times Ionized Sodium (NaIX)")
        ENUM_VAL(Species, NineTimesIonizedSodium, "Nine times Ionized Sodium (NaX)")
        ENUM_VAL(Species, TenTimesIonizedSodium, "Ten times Ionized Sodium (NaXI)")

        ENUM_VAL(Species, IonizedMagnesium, "Ionized Magnesium (MgII)")
        ENUM_VAL(Species, DoublyIonizedMagnesium, "Doubly Ionized Magnesium (MgIII)")
        ENUM_VAL(Species, TriplyIonizedMagnesium, "Triply Ionized Magnesium (MgIV)")
        ENUM_VAL(Species, FourTimesIonizedMagnesium, "Four times Ionized Magnesium (MgV)")
        ENUM_VAL(Species, FiveTimesIonizedMagnesium, "Five times Ionized Magnesium (MgVI)")
        ENUM_VAL(Species, SixTimesIonizedMagnesium, "Six times Ionized Magnesium (MgVII)")
        ENUM_VAL(Species, SevenTimesIonizedMagnesium, "Seven times Ionized Magnesium (MgVIII)")
        ENUM_VAL(Species, EightTimesIonizedMagnesium, "Eight times Ionized Magnesium (MgIX)")
        ENUM_VAL(Species, NineTimesIonizedMagnesium, "Nine times Ionized Magnesium (MgX)")
        ENUM_VAL(Species, TenTimesIonizedMagnesium, "Ten times Ionized Magnesium (MgXI)")

        ENUM_VAL(Species, AtomicSilicon, "Atomic Silicon (SiI)")
        ENUM_VAL(Species, IonizedSilicon, "Ionized Silicon (SiII)")
        ENUM_VAL(Species, DoublyIonizedSilicon, "Doubly Ionized Silicon (SiIII)")
        ENUM_VAL(Species, TriplyIonizedSilicon, "Triply Ionized Silicon (SiIV)")
        ENUM_VAL(Species, FourTimesIonizedSilicon, "Four times Ionized Silicon (SiV)")
        ENUM_VAL(Species, FiveTimesIonizedSilicon, "Five times Ionized Silicon (SiVI)")
        ENUM_VAL(Species, SixTimesIonizedSilicon, "Six times Ionized Silicon (SiVII)")
        ENUM_VAL(Species, SevenTimesIonizedSilicon, "Seven times Ionized Silicon (SiVIII)")
        ENUM_VAL(Species, EightTimesIonizedSilicon, "Eight times Ionized Silicon (SiIX)")
        ENUM_VAL(Species, NineTimesIonizedSilicon, "Nine times Ionized Silicon (SiX)")
        ENUM_VAL(Species, TenTimesIonizedSilicon, "Ten times Ionized Silicon (SiXI)")

        ENUM_VAL(Species, AtomicSulfur, "Atomic Sulfur (SI)")
        ENUM_VAL(Species, IonizedSulfur, "Ionized Sulfur (SII)")
        ENUM_VAL(Species, DoublyIonizedSulfur, "Doubly Ionized Sulfur (SIII)")
        ENUM_VAL(Species, TriplyIonizedSulfur, "Triply Ionized Sulfur (SIV)")
        ENUM_VAL(Species, FourTimesIonizedSulfur, "Four times Ionized Sulfur (SV)")
        ENUM_VAL(Species, FiveTimesIonizedSulfur, "Five times Ionized Sulfur (SVI)")
        ENUM_VAL(Species, SixTimesIonizedSulfur, "Six times Ionized Sulfur (SVII)")
        ENUM_VAL(Species, SevenTimesIonizedSulfur, "Seven times Ionized Sulfur (SVIII)")
        ENUM_VAL(Species, EightTimesIonizedSulfur, "Eight times Ionized Sulfur (SIX)")
        ENUM_VAL(Species, NineTimesIonizedSulfur, "Nine times Ionized Sulfur (SIX)")
        ENUM_VAL(Species, TenTimesIonizedSulfur, "Ten times Ionized Sulfur (SXI)")

        ENUM_VAL(Species, DoublyIonizedArgon, "Doubly Ionized Argon (ArIII)")
        ENUM_VAL(Species, TriplyIonizedArgon, "Triply Ionized Argon (ArIV)")
        ENUM_VAL(Species, FourTimesIonizedArgon, "Four times Ionized Argon (ArV)")
        ENUM_VAL(Species, FiveTimesIonizedArgon, "Five times Ionized Argon (ArVI)")
        ENUM_VAL(Species, SixTimesIonizedArgon, "Six times Ionized Argon (ArVII)")
        ENUM_VAL(Species, SevenTimesIonizedArgon, "Seven times Ionized Argon (ArVIII)")
        ENUM_VAL(Species, EightTimesIonizedArgon, "Eight times Ionized Argon (ArIX)")
        ENUM_VAL(Species, NineTimesIonizedArgon, "Nine times Ionized Argon (ArX)")
        ENUM_VAL(Species, TenTimesIonizedArgon, "Ten times Ionized Argon (ArXI)")

        ENUM_VAL(Species, IonizedIron, "Ionized Iron (FeII)")
        ENUM_VAL(Species, DoublyIonizedIron, "Doubly Ionized Iron (FeIII)")
        ENUM_VAL(Species, TriplyIonizedIron, "Triply Ionized Iron (FeIV)")
        ENUM_VAL(Species, FourTimesIonizedIron, "Four times Ionized Iron (FeV)")
        ENUM_VAL(Species, FiveTimesIonizedIron, "Five times Ionized Iron (FeVI)")
        ENUM_VAL(Species, SixTimesIonizedIron, "Six times Ionized Iron (FeVII)")
        ENUM_VAL(Species, SevenTimesIonizedIron, "Seven times Ionized Iron (FeVIII)")
        ENUM_VAL(Species, EightTimesIonizedIron, "Eight times Ionized Iron (FeIX)")
        ENUM_VAL(Species, NineTimesIonizedIron, "Nine times Ionized Iron (FeX)")
        ENUM_VAL(Species, TenTimesIonizedIron, "Ten times Ionized Iron (FeXI)")

    ENUM_END()

    /** The enumeration type indicating the assumed initial level population. */
    ENUM_DEF(InitialLevelPopsCase, LTE, CollisionallyExcited, Custom)
        ENUM_VAL(InitialLevelPopsCase, LTE, "LTE (the detailed balance condition between collisional transitions)")
        ENUM_VAL(InitialLevelPopsCase, CollisionallyExcited,
                 "the detailed balance condition between collisional transitions and spontaneous emission")
        ENUM_VAL(InitialLevelPopsCase, Custom, "custom level populations provided in a file")
    ENUM_END()

    ITEM_CONCRETE(NonLTELineGasMix, EmittingGasMix,
                  "A gas mix supporting rotational transitions in specific molecules and atoms")
        ATTRIBUTE_TYPE_INSERT(NonLTELineGasMix, "CustomMediumState,DynamicState")

        PROPERTY_ENUM(species, Species, "the molecular or atomic species being represented")
        ATTRIBUTE_DEFAULT_VALUE(species, "CarbonMonoxide")

        PROPERTY_INT(numEnergyLevels, "the number of energy levels used (or 999 for all supported)")
        ATTRIBUTE_MIN_VALUE(numEnergyLevels, "2")
        ATTRIBUTE_MAX_VALUE(numEnergyLevels, "999")
        ATTRIBUTE_DEFAULT_VALUE(numEnergyLevels, "999")
        ATTRIBUTE_DISPLAYED_IF(numEnergyLevels, "Level2")

        PROPERTY_DOUBLE(defaultTemperature, "the default temperature of the gas")
        ATTRIBUTE_QUANTITY(defaultTemperature, "temperature")
        ATTRIBUTE_MIN_VALUE(defaultTemperature, "[0")
        ATTRIBUTE_MAX_VALUE(defaultTemperature, "1e9]")
        ATTRIBUTE_DEFAULT_VALUE(defaultTemperature, "1000")
        ATTRIBUTE_DISPLAYED_IF(defaultTemperature, "Level2")

        PROPERTY_DOUBLE_LIST(defaultCollisionPartnerRatios,
                             "the default relative abundances of the collisional partners")
        ATTRIBUTE_MIN_VALUE(defaultCollisionPartnerRatios, "[0")
        ATTRIBUTE_MAX_VALUE(defaultCollisionPartnerRatios, "1e50]")
        ATTRIBUTE_DEFAULT_VALUE(defaultCollisionPartnerRatios, "1e4")
        ATTRIBUTE_DISPLAYED_IF(defaultCollisionPartnerRatios, "Level2")

        PROPERTY_DOUBLE(defaultTurbulenceVelocity, "the default (non-thermal) turbulence velocity")
        ATTRIBUTE_QUANTITY(defaultTurbulenceVelocity, "velocity")
        ATTRIBUTE_MIN_VALUE(defaultTurbulenceVelocity, "[0 km/s")
        ATTRIBUTE_MAX_VALUE(defaultTurbulenceVelocity, "100000 km/s]")
        ATTRIBUTE_DEFAULT_VALUE(defaultTurbulenceVelocity, "0 km/s")
        ATTRIBUTE_DISPLAYED_IF(defaultTurbulenceVelocity, "Level2")

        PROPERTY_DOUBLE(maxChangeInLevelPopulations,
                        "the maximum relative change for the level populations in a cell to be considered converged")
        ATTRIBUTE_MIN_VALUE(maxChangeInLevelPopulations, "[0")
        ATTRIBUTE_MAX_VALUE(maxChangeInLevelPopulations, "1]")
        ATTRIBUTE_DEFAULT_VALUE(maxChangeInLevelPopulations, "0.05")
        ATTRIBUTE_DISPLAYED_IF(maxChangeInLevelPopulations, "Level2")

        PROPERTY_DOUBLE(maxFractionNotConvergedCells,
                        "the maximum fraction of not-converged cells for all cells to be considered converged")
        ATTRIBUTE_MIN_VALUE(maxFractionNotConvergedCells, "[0")
        ATTRIBUTE_MAX_VALUE(maxFractionNotConvergedCells, "1]")
        ATTRIBUTE_DEFAULT_VALUE(maxFractionNotConvergedCells, "0.01")
        ATTRIBUTE_DISPLAYED_IF(maxFractionNotConvergedCells, "Level2")

        PROPERTY_DOUBLE(maxChangeInGlobalLevelPopulations,
                        "the maximum relative change for the global level populations to be considered converged")
        ATTRIBUTE_MIN_VALUE(maxChangeInGlobalLevelPopulations, "[0")
        ATTRIBUTE_MAX_VALUE(maxChangeInGlobalLevelPopulations, "1]")
        ATTRIBUTE_DEFAULT_VALUE(maxChangeInGlobalLevelPopulations, "0.05")
        ATTRIBUTE_DISPLAYED_IF(maxChangeInGlobalLevelPopulations, "Level2")

        PROPERTY_ENUM(initialLevelPopsCase, InitialLevelPopsCase, "the assumed initial level population")
        ATTRIBUTE_DEFAULT_VALUE(initialLevelPopsCase, "LTE")
        ATTRIBUTE_DISPLAYED_IF(initialLevelPopsCase, "Level2")

        PROPERTY_DOUBLE(lowestOpticalDepth, "Lower limit of (negative) optical depth along a cell diagonal")
        ATTRIBUTE_MIN_VALUE(lowestOpticalDepth, "[-10")
        ATTRIBUTE_MAX_VALUE(lowestOpticalDepth, "0]")
        ATTRIBUTE_DEFAULT_VALUE(lowestOpticalDepth, "-2")
        ATTRIBUTE_DISPLAYED_IF(lowestOpticalDepth, "Level3")

        PROPERTY_BOOL(storeMeanIntensities, "store the mean radiation field intensity at each transition line")
        ATTRIBUTE_DEFAULT_VALUE(storeMeanIntensities, "false")
        ATTRIBUTE_DISPLAYED_IF(storeMeanIntensities, "Level3")

        PROPERTY_STRING(initialLevelPopsFilename, "the name of the file with initial level populations")
        ATTRIBUTE_REQUIRED_IF(initialLevelPopsFilename, "false")
        ATTRIBUTE_DISPLAYED_IF(initialLevelPopsFilename, "Level3")

        PROPERTY_BOOL(updateDynamicStatesFlag, "update level populations from the initial one")
        ATTRIBUTE_DEFAULT_VALUE(updateDynamicStatesFlag, "true")
        ATTRIBUTE_DISPLAYED_IF(updateDynamicStatesFlag, "Level3")

        PROPERTY_BOOL(errorForGaussianIntegral,
                  "treat insufficient Gaussian-integral accuracy as a fatal error (otherwise issue a warning)")
        ATTRIBUTE_DEFAULT_VALUE(errorForGaussianIntegral, "true")
        ATTRIBUTE_DISPLAYED_IF(errorForGaussianIntegral, "Level3")

        PROPERTY_DOUBLE(lowestBranchingRatio, "the minimum branching ratio for a transition to be considered")
        ATTRIBUTE_MIN_VALUE(lowestBranchingRatio, "[0")
        ATTRIBUTE_MAX_VALUE(lowestBranchingRatio, "1]")
        ATTRIBUTE_DEFAULT_VALUE(lowestBranchingRatio, "0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function loads the required information on the configured species from resource files:

        - the molecular weight of the species;

        - the supported energy levels with their corresponding degeneracies;

        - the Einstein A coefficients for the radiative transitions between the energy levels;

        - the gas temperature grid points used for discretizing the collisional coefficients K;

        - the excitation coefficients K for the collisional transitions between the energy levels,
        given for each of the temperature grid points.

        The last two items are repeated for each collisional interaction partner of the species.

        In a second step, the function calculates and caches some further values, including for
        example the Einstein B coefficients. */
    void setupSelfBefore() override;

    //======== Capabilities =======

public:
    /** This function returns true, indicating that this material may have a negative absorption
        cross section. This happens when inverted level populations cause net stimulated emission.
        */
    bool hasNegativeExtinction() const override;

    /** This function returns true, indicating that the cross sections returned by this material
        mix depend on the values of specific state variables other than the number density. */
    bool hasExtraSpecificState() const override;

    /** This function returns DynamicStateType::PrimaryIfMergedIterations, indicating that this
        material mix has a dynamic medium state with updates that are considered to affect primary
        emission when the simulation has merged iterations, and only secondary emission if not. */
    DynamicStateType hasDynamicMediumState() const override;

    /** This function returns true, indicating that this material supports secondary line emission
        from gas. */
    bool hasLineEmission() const override;

    //======== Medium state setup =======

public:
    /** This function returns the number and type of import parameters required by this particular
        material mix as a list of SnapshotParameter objects. For this class, the function returns a
        descriptor for the number densities of the collisional partners and for the turbulence
        velocity. Importing the kinetic gas temperature should be enabled through the corresponding
        standard configuration flag. */
    vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns a list of StateVariable objects describing the specific state
        variables used by the receiving material mix. For this class, the function returns a list
        containing descriptors for the properties defined in the input model and for a number of
        variables to hold the level populations derived from the radiation field and related
        information. */
    vector<StateVariable> specificStateVariableInfo() const override;

    /** This function initializes the specific state variables requested by this material mix
        through the specificStateVariableInfo() function except for the number density of the
        species under consideration. For this class, the function uses the imported values, or if
        not available, the user-configured default values. The level populations are left at zero.
        */
    void initializeSpecificState(MaterialState* state, double metallicity, double temperature,
                                 const Array& params) const override;

    //======== Medium state updates =======

    /** Based on the specified radiation field and the input model properties found in the given
        material state, this function determines the level populations for the supported
        transitions and stores these results back in the given material state. The function returns
        the update status as described for the UpdateStatus class. */
    UpdateStatus updateSpecificState(MaterialState* state, const Array& Jv) const override;

    /** This function returns true if the state of the medium component corresponding to this
        material mix can be considered to be converged based on the given spatial cell statistics
        and aggregate material states, and false otherwise. For more information on the convergence
        criteria, see the corresponding section in the class header. */
    bool isSpecificStateConverged(int numCells, int numUpdated, int numNotConverged, MaterialState* currentAggregate,
                                  MaterialState* previousAggregate) const override;

    //======== Low-level material properties =======

public:
    /** This function returns the mass of a molecule or atom of the species under consideration. */
    double mass() const override;

    /** This function should return the absorption cross section using default properties. Because
        this value is hard to calculate for this material mix, this function returns zero. */
    double sectionAbs(double lambda) const override;

    /** This function should return the scattering cross section using default properties. Because
        this value is hard to calculate for this material mix, this function returns zero. */
    double sectionSca(double lambda) const override;

    /** This function should return the extinction cross section using default properties. Because
        this value is hard to calculate for this material mix, this function returns zero. */
    double sectionExt(double lambda) const override;

    //======== High-level photon life cycle =======

    /** This function returns the absorption opacity \f$k^\text{abs}= n_\mathrm{mol}
        \varsigma^\text{abs}\f$ for the given wavelength and material state. The photon packet
        properties are not used. */
    double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the scattering opacity \f$k^\text{sca}\f$ which is trivially zero at
        all wavelengths. */
    double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the extinction opacity \f$k^\text{ext}=k^\text{abs}\f$ for the given
        wavelength and material state, which equals the absorption opacity because the scattering
        opacity is zero. The photon packet properties are not used. */
    double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function does nothing because the lines under consideration do not scatter. */
    bool peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs, Direction bfky,
                           const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function does nothing because the lines under consideration do not scatter. */
    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======== Secondary emission =======

    /** This function returns a list with the line centers of the supported transitions. */
    Array lineEmissionCenters() const override;

    /** This function returns a list with the masses of the particles emitting each of the lines,
        i.e. the mass of a molecule or atom of the species under consideration in each case. */
    Array lineEmissionMasses() const override;

    /** This function returns a list with the line luminosities for the supported transitions in
        the spatial cell and medium component represented by the specified material state and the
        receiving material mix when it would be embedded in the specified radiation field. */
    Array lineEmissionSpectrum(const MaterialState* state, const Array& Jv) const override;

    //======== Temperature =======

    /** This function returns an indicative temperature of the material mix when it would be
        embedded in a given radiation field. The implementation in this class ignores the radiation
        field and returns the temperature stored in the specific state for the relevant spatial
        cell and medium component. This value corresponds to the temperature defined by the input
        model at the start of the simulation. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //======================== Data Members ========================

private:
    // Data members loaded from text resource files and/or precalculated in setupSelfBefore()
    // (only the energy levels and transitions actually used are stored in the data members)

    string _name;  // human readable species name

    // mass
    double _mass{0.};  // particle mass for the species under consideration

    // energy levels
    int _numLevels{0};       // the number of energy levels -- index p
    vector<double> _energy;  // the energy of each energy level
    vector<double> _weight;  // the weight (degeneracy) of each energy level

    // radiative transitions
    int _numLines{0};             // the number of radiative transitions -- index k
    vector<int> _indexUpRad;      // the upper energy level index for each radiative transition
    vector<int> _indexLowRad;     // the lower energy level index for each radiative transition
    vector<double> _einsteinA;    // the Einstein A coefficient for each radiative transition
    vector<double> _einsteinBul;  // the Einstein Bul coefficient for each radiative transition
    vector<double> _einsteinBlu;  // the Einstein Blu coefficient for each radiative transition
    vector<double> _branchRatio;  // the branching ratio for each radiative transition from the same upper energy level
    Array _center;                // the central emission wavelength for each radiative transition

    // collisional transitions
    struct ColPartner  // data structure holding information on a collisional partner
    {
        string name;              // human readable species name
        Array T;                  // the temperature grid points
        int numColTrans{0};       // the number of collisional transitions -- index t
        vector<int> indexUpCol;   // the upper energy level index for each collisional transition
        vector<int> indexLowCol;  // the lower energy level index for each collisional transition
        vector<Array> Kul;        // the coefficient for each collisional transition and for each temperature
    };
    int _numColPartners{0};          // the number of collisional interaction partners -- index c
    vector<ColPartner> _colPartner;  // the data for each collisional partner

    // the radiation field wavelength grid for this simulation
    int _numWavelengths{0};  // the number of wavelength bins -- index ell
    Array _lambdav;          // characteristic wavelengths
    Array _dlambdav;         // wavelength bin widths

    // imported initial level populations (technical expert option)
    vector<Array> _initLevelPops;  // initial level populations for each cell -- indices m, p

private:
    // Data members indicating custom variable indices; initialized in specificStateVariableInfo()
    int _indexKineticTemperature{0};
    int _indexFirstColPartnerDensity{0};
    int _indexFirstLevelPopulation{0};
    int _indexFirstMeanIntensity{0};
};

////////////////////////////////////////////////////////////////////

#endif
