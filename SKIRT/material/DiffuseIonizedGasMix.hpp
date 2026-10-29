/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DIFFUSEIONIZEDGASMIX_HPP
#define DIFFUSEIONIZEDGASMIX_HPP

#include "EmittingGasMix.hpp"
#include "GasLineEmission.hpp"
#include "PhotoIonizationSolver.hpp"
#include "StoredTable.hpp"
#include <array>

//////////////////////////////////////////////////////////////////////

/** An instance of the DiffuseIonizedGasMix class represents diffuse ionized gas regions (HII regions,
    diffuse ionized gas, ionization fronts) for radiative transfer simulations. The module uses a
    hybrid approach that combines pre-computed Cloudy tables with an inline photoionization solver.

    <b>Applicability and naming note</b>

    The class name "DiffuseIonizedGasMix" reflects the module's original motivating application:
    modelling the diffuse ionised medium in cosmological galaxy simulations, where dense HII
    regions are often unresolved and handled by sub-grid prescriptions (Kapoor et al. 2026a; the
    multi-metallicity and per-cell-abundance extensions are presented in Kapoor et al. 2026b).
    The module itself is general-purpose: it operates on the local radiation field and gas
    density in each cell, agnostic to the specific astrophysical context.

    The STAB tables extend to n_H = 10^5 cm^-3 and log U from -6.5 to +1.5, making the module
    applicable across the photoionised-plasma regime, from the diffuse ionised medium
    (n_H << 1 cm^-3) through individual HII regions, planetary nebulae, and dense star-forming
    clumps (n_H up to 10^5 cm^-3).

    The STAB tables for temperature and opacity are pre-computed by running Cloudy on a grid of
    5-bin step-function SEDs that are zero above 6 Ryd; the lookup parameters are the
    ionisation parameter log U (integrated from 1 Ryd over the full simulation wavelength grid)
    and four spectral-shape ratios R2-R5 covering 1.0-6.0 Ryd. The inline ionisation balance,
    by contrast, is solved by PhotoIonizationSolver directly from the full radiation field
    using Verner photoionisation cross-sections extending to keV. As a result, the ion balance
    accommodates any incident SED, but the STAB-based T and opacity ignore any flux above
    6 Ryd. When the radiation-field wavelength grid extends beyond 6 Ryd and the incident SED
    carries significant flux there (e.g. AGN power-law spectra into X-rays), the T and opacity
    lookups become approximate; the size of the bias scales with the fraction of ionising flux
    that sits in that regime.

    - <b>Temperature and opacity</b> are determined from pre-computed STAB tables based on
      the ionization parameter (log U), a 5-bin spectral characterization (R2-R5),
      metallicity, and gas density.
    - <b>Ionization balance</b> is solved every iteration by the inline PhotoIonizationSolver
      at the STAB-derived temperature, yielding ion fractions for all tracked species (H, He,
      C, N, O, Ne, Mg, Si, S, Fe) and the electron density. Key ion fractions (x_HII, x_NII,
      x_OI, x_OII, x_OIII, x_SII) and n_e are stored as custom state variables for diagnostics.
    - <b>Emission</b> is computed from the converged radiation field, temperature, and ion
      fractions: GasLineEmission (this class) computes the line luminosities, while
      GasContinuumEmission computes the nebular continuum (free-bound, free-free, two-photon). By
      default (\em includeExtendedLines=true), H and He recombination lines from Case B are obtained
      from emissivity tables, and collisional lines from the statistical-equilibrium solver on the
      atomic data. This obviously requires the corresponding resources to be installed. If \em
      includeExtendedLines is set to false, the Storey & Hummer 1995 P_B form for H and precomputed
      CHIANTI q_col tables are used as a fallback.
    - <b>Diffuse reemission</b> of absorbed ionizing photons follows Wood, Mathis & Ercolano (2004),
      with pre-tabulated temperature-dependent spectra for H/He Lyman continuum and He two-photon.

    <b>5-Bin Radiation Field Characterization</b>

    The ionizing radiation field is characterized using 5 energy bins:
    - Bin 1: 1.0 - 1.8 Ryd (912 - 506 A)
    - Bin 2: 1.8 - 2.58 Ryd (506 - 353 A)
    - Bin 3: 2.58 - 3.52 Ryd (353 - 259 A)
    - Bin 4: 3.52 - 4.0 Ryd (259 - 228 A)
    - Bin 5: 4.0 - 6.0 Ryd (228 - 152 A)

    The spectral shape is described by four ratios relative to the first bin:
    R2 = J2/J1, R3 = J3/J1, R4 = J4/J1, R5 = J5/J1

    <b>Dual-Table System</b>

    Separate STAB tables for temperature and opacity optimize accuracy across different
    ionization regimes:

    - Standard Table: log U sampling (-4.0, -3.3, -2.6, -2.0, -1.3, -0.6, 0.1, 0.8, 1.5)
      for general ISM (log U >= -4)
    - Transition Table: log U sampling (-6.5, -6.1, -5.7, -5.2, -4.8, -4.4, -4.0)
      for ionization fronts (log U <= -4), capturing radiation hardening at nebular edges.
      The transition table uses broader shape-ratio ranges to capture the harder spectra
      near ionization fronts.

    The appropriate table is automatically selected based on log U, with smooth blending in
    overlap regions.

    <b>Required Input</b>

    The input model must provide spatial distributions of:
    - Number density
    - Metallicity

    <b>STAB Table Structure</b>

    - Temperature: 7D (Z, n_H, logU, log10(R2), log10(R3), log10(R4), log10(R5))
    - Opacity: 9D (wavelength, logU, log10(R2), log10(R3), log10(R4), log10(R5), Z, n_H, index)
    - Reemission spectra: 3D (spectrumType, T, wavelength) -> cumulative probability

    <b>Convergence</b>

    Convergence of the primary emission iterations is governed by three criteria. The per-cell
    criterion OR the plateau criterion must pass, AND the global criterion must pass:
    - Per-cell: each cell's relative change in temperature |delta T|/T and ionization parameter
      |delta logU|/|logU| are compared to maxChangeInIonization. A cell is "not converged" if
      either exceeds the threshold. The fraction of not-converged cells must be below
      maxFractionNotConvergedCells.
    - Plateau: the converged cell fraction is tracked over the last 3 iterations. If the change
      between consecutive values is below stabilityConvergenceThreshold (a fraction, displayed as
      percent in the log), convergence is declared even if the per-cell criterion is not met. This
      handles cases where Monte Carlo noise prevents the per-cell fraction from reaching the target.
    - Global: the relative change |delta sum(nHp*V)| / sum(nHp*V) in total ionized hydrogen must
      be below maxChangeInGlobalIonizedH. This prevents false convergence when trivially-neutral
      cells dominate the per-cell count while the ionization structure is still evolving.

    All threshold properties are specified as fractions (0-1); log output displays percentages.

    <b>Gas density</b>

    The number density of this mix is the number density of hydrogen plus helium nuclei,
    n = n_H (1 + y_He), with y_He = n_He/n_H the helium abundance of the cell (the imported
    column in PerCell mode, the value implied by the metallicity in SolarScaled mode). The mix
    derives the hydrogen density as n_H = n / (1 + y_He), so a snapshot that lists n_H must be
    multiplied by (1 + y_He) before it is imported as a number density. Accordingly, the mass per
    particle is the mean mass per nucleus, m_p (1 + 4 y_He) / (1 + y_He). Because it is a property
    of the mix rather than of a cell, it uses the y_He of the \em defaultMetallicity value; it
    enters only conversions between number and mass density (for example reported gas masses, or
    a snapshot that lists mass density), not the gas physics.

    <b>Density Ceiling</b>

    If maxHydrogenDensity is set to a positive value, the hydrogen number density n_H is clamped
    to that ceiling (in cm^-3) at the start of each cell update, before it is used anywhere:
    ionization parameter (logU), ionization fractions, temperature/opacity STAB lookups, and
    emission. This ensures self-consistency across all calculations. The default value is
    10^5 cm^-3, matching the upper edge of the stab density axis; setting it to 0 disables
    the ceiling.

    For opacity, cells below the minimum STAB table density are linearly scaled
    (opacity proportional to n). Cells above the maximum table density use StoredTable clamping.
    Cells with log U at or below -6.5, the lower edge of the tables, use the Verner bound-free
    opacity computed from the ion fractions of the ionization balance instead of the table, as all
    cells do when \em useCloudyOpacity is false.

    <b>Radiation field wavelength grid</b>

    The field averages over the five energy bins (R2-R5) and the ionizing photon flux (log U) are
    computed from the mean intensity stored per radiation field bin, and are exact when the grid has
    bin borders at the energy bin boundaries (1, 1.8, 2.58, 3.52, 4 and 6 Ryd). The ionization
    balance evaluates each bin at its characteristic wavelength, so borders at the ionization edges
    of the tracked ions, and in particular at the hydrogen edge, avoid counting a bin as wholly
    above or below an edge. During setup, the mix warns when an energy bin boundary or the hydrogen
    edge falls inside a bin, and lists the other ionization edges that do.

    Two radiation field wavelength grids are recommended for this mix. Both have borders at the six
    energy bin boundaries and at ionization edges of the tracked ions, as specified below, and
    enough bins in between that the photoionization rates, the hydrogen and helium photoheating and
    the absorption in the transport are within 2 %, and the power absorbed by dust within 1 %, of a
    fine reference grid, for blackbody fields of 30 to 50 kK behind H I columns up to
    \f$10^{18}\,\mathrm{cm}^{-2}\f$. The medium grid (70 bins) has borders at the edges of H, He, C,
    N, O, Ne and S and is the default choice. The fine grid (83 bins) adds the edges of Mg, Si and
    Fe, so that every ionization edge of the tracked ions that falls within the grid has a border;
    use it when the lines of these elements matter. Both grids start at 0.0092 micron, so edges
    above 135 eV lie outside them and have no border. Save the borders of the chosen grid as a text
    file and use it as:

    \verbatim
    <radiationFieldWLG type="DisjointWavelengthGrid">
        <FileBorderWavelengthGrid filename="dig_rf_grid_medium.txt" characteristic="Logarithmic"/>
    </radiationFieldWLG>
    \endverbatim

    The borders match the band boundaries and edge energies used in this class; the setup check
    reports any border that no longer matches them.

    Medium grid (71 borders, micron):

    \verbatim
    0.0091882500
    0.0098244214
    0.0103412919
    0.0106098368
    0.0108853554
    0.0117418171
    0.0121949959
    0.0126656653
    0.0127660830
    0.0139243820
    0.0151877764
    0.0155906237
    0.0160041562
    0.0160165609
    0.0175477571
    0.0183674044
    0.0192253370
    0.0195373776
    0.0202543879
    0.0209977120
    0.0217683158
    0.0221641589
    0.0225672002
    0.0227816647
    0.0242853361
    0.0250739902
    0.0258882553
    0.0261294411
    0.0262067635
    0.0302695796
    0.0353030177
    0.0353204103
    0.0355969561
    0.0386139068
    0.0418865534
    0.0438740698
    0.0459558938
    0.0481365004
    0.0492653033
    0.0498395927
    0.0504205767
    0.0506259215
    0.0508548802
    0.0531436768
    0.0575065853
    0.0634692500
    0.0700501632
    0.0773134293
    0.0853297990
    0.0910309824
    0.0911266587
    0.0911648515
    0.1101103002
    0.1196758668
    0.1587769964
    0.1828849982
    0.2106534530
    0.2426381481
    0.2794792589
    0.3707921945
    0.4919393735
    0.6526683970
    1.1488268704
    1.5241776371
    2.0221649835
    3.5594146741
    6.2652814808
    7.2165743031
    8.3123072493
    9.5744114736
    11.0281480600
    \endverbatim

    Fine grid (84 borders, micron):

    \verbatim
    0.0091882500
    0.0098244214
    0.0099187358
    0.0103908116
    0.0106352093
    0.0108853554
    0.0113434765
    0.0119153580
    0.0125160709
    0.0126656653
    0.0127660830
    0.0139243820
    0.0151877764
    0.0154709506
    0.0160041562
    0.0160165609
    0.0165290225
    0.0178262735
    0.0185125934
    0.0192253370
    0.0195373776
    0.0202543879
    0.0209977120
    0.0217683158
    0.0221641589
    0.0225672002
    0.0226248536
    0.0227816647
    0.0242853361
    0.0250739902
    0.0258882553
    0.0261294411
    0.0262067635
    0.0274665924
    0.0302695796
    0.0353030177
    0.0353204103
    0.0355969561
    0.0370212595
    0.0386984459
    0.0404516144
    0.0418865534
    0.0459558938
    0.0481365004
    0.0492653033
    0.0498395927
    0.0504205767
    0.0506259215
    0.0508548802
    0.0531436768
    0.0575065853
    0.0660363530
    0.0707645632
    0.0758313138
    0.0765807276
    0.0824363019
    0.0853297990
    0.0910309824
    0.0911266587
    0.0911648515
    0.1101103002
    0.1196758668
    0.1349131784
    0.1432446036
    0.1520905275
    0.1569023007
    0.1621556343
    0.1850119894
    0.2110900207
    0.2408438339
    0.2747915422
    0.3577165393
    0.4656661608
    0.6061921927
    0.7891253550
    1.3372645004
    1.7408164215
    2.2661498997
    3.8402540158
    6.5077561320
    7.4250451627
    8.4716290147
    9.6657322063
    11.0281480600
    \endverbatim

    <b>Abundances and gas-phase depletion</b>

    The metallicity and per-cell element abundances supplied to this mix are
    interpreted as total values, including atoms that are locked up in dust.
    The temperature and opacity tables shipped with this mix were generated by
    Cloudy runs in which a fraction of the heavy elements was placed in the dust
    phase following the recipe of Gunasekera et al. (2023), based on the
    Jenkins (2009) depletion pattern at F* = 0.5 with no depletion applied to
    sulphur. The line-emission and ionization-balance step inside
    this mix applies that same gas-phase fraction to the abundance vector before
    it is used, so the analytical step is consistent with the physics in the
    tables. Inputs must therefore be total abundances; subtracting dust before
    passing them in would apply the depletion twice.

    The abundanceMode enum selects how metal and helium abundances are handed to
    the analytical PIO/emission solver: SolarScaled uses the same prescription
    the stabs were built with (solar pattern * Zfrac * g(dN,dC) for metals,
    Gutkin Y(Z) for helium), so emission is fully self-consistent with the
    tabulated T and kappa. PerCell uses the snapshot's abundances in the
    emission calculation via 9 imported columns (8 metals + He); the stab
    T/kappa stay on the Gutkin baseline.
*/
class DiffuseIonizedGasMix : public EmittingGasMix
{
    /** Selects how metal and helium abundances are handed to the analytical
        emission solver. See the class documentation section "Abundances and
        gas-phase depletion" for the full description. */
    ENUM_DEF(AbundanceMode, SolarScaled, PerCell)
        ENUM_VAL(AbundanceMode, SolarScaled, "solar pattern scaled by metallicity")
        ENUM_VAL(AbundanceMode, PerCell, "per-cell abundances imported from the snapshot")
    ENUM_END()

    ITEM_CONCRETE(DiffuseIonizedGasMix, EmittingGasMix, "a diffuse ionized gas mix including emission")
        ATTRIBUTE_TYPE_INSERT(DiffuseIonizedGasMix, "CustomMediumState,DynamicState")

        PROPERTY_DOUBLE(defaultMetallicity, "the default metallicity of the gas")
        ATTRIBUTE_MIN_VALUE(defaultMetallicity, "]0")
        ATTRIBUTE_MAX_VALUE(defaultMetallicity, "0.2]")
        ATTRIBUTE_DEFAULT_VALUE(defaultMetallicity, "0.02")
        ATTRIBUTE_DISPLAYED_IF(defaultMetallicity, "Level2")

        PROPERTY_DOUBLE(defaultTemperature, "the fixed temperature of the gas")
        ATTRIBUTE_QUANTITY(defaultTemperature, "temperature")
        ATTRIBUTE_MIN_VALUE(defaultTemperature, "[1000")
        ATTRIBUTE_MAX_VALUE(defaultTemperature, "20000]")
        ATTRIBUTE_DEFAULT_VALUE(defaultTemperature, "10000")
        ATTRIBUTE_DISPLAYED_IF(defaultTemperature, "Level2")

        PROPERTY_BOOL(useCloudyTemperature,
                      "use the Cloudy-tabulated temperature (disabling keeps T fixed at defaultTemperature)")
        ATTRIBUTE_DEFAULT_VALUE(useCloudyTemperature, "true")
        ATTRIBUTE_DISPLAYED_IF(useCloudyTemperature, "Level3")

        PROPERTY_BOOL(useCloudyOpacity,
                      "use the Cloudy-tabulated opacity (disabling falls back to Verner bound-free opacity)")
        ATTRIBUTE_DEFAULT_VALUE(useCloudyOpacity, "true")
        ATTRIBUTE_DISPLAYED_IF(useCloudyOpacity, "Level3")

        // ===== Abundance specification (metals + helium, for analytical emission solver) =====
        // T/kappa stabs always use the Gutkin baseline regardless of this setting.
        PROPERTY_ENUM(abundanceMode, AbundanceMode, "abundance specification for the emission solver")
        ATTRIBUTE_DEFAULT_VALUE(abundanceMode, "SolarScaled")
        ATTRIBUTE_DISPLAYED_IF(abundanceMode, "Level2")

        // Range matches the (delta_N, delta_C) cross sampled in the multiZ stab DeltaMap:
        // delta_N in [-1.5, +1.0], delta_C in [-1.0, +1.0]. Values outside the cross are
        // endpoint-clamped at runtime.
        PROPERTY_DOUBLE(solarScalingN, "nitrogen scaling factor on solar pattern * Zfrac (SolarScaled mode)")
        ATTRIBUTE_MIN_VALUE(solarScalingN, "[0.03162")
        ATTRIBUTE_MAX_VALUE(solarScalingN, "10.0]")
        ATTRIBUTE_DEFAULT_VALUE(solarScalingN, "1")
        ATTRIBUTE_RELEVANT_IF(solarScalingN, "abundanceModeSolarScaled")
        ATTRIBUTE_DISPLAYED_IF(solarScalingN, "Level2")

        PROPERTY_DOUBLE(solarScalingC, "carbon scaling factor on solar pattern * Zfrac (SolarScaled mode)")
        ATTRIBUTE_MIN_VALUE(solarScalingC, "[0.1")
        ATTRIBUTE_MAX_VALUE(solarScalingC, "10.0]")
        ATTRIBUTE_DEFAULT_VALUE(solarScalingC, "1")
        ATTRIBUTE_RELEVANT_IF(solarScalingC, "abundanceModeSolarScaled")
        ATTRIBUTE_DISPLAYED_IF(solarScalingC, "Level2")

        PROPERTY_DOUBLE(maxChangeInIonization, "maximum relative change in T or logU per cell for convergence")
        ATTRIBUTE_MIN_VALUE(maxChangeInIonization, "[0")
        ATTRIBUTE_MAX_VALUE(maxChangeInIonization, "1]")
        ATTRIBUTE_DEFAULT_VALUE(maxChangeInIonization, "0.01")
        ATTRIBUTE_DISPLAYED_IF(maxChangeInIonization, "Level2")

        PROPERTY_DOUBLE(maxFractionNotConvergedCells, "maximum fraction of cells allowed to not converge")
        ATTRIBUTE_MIN_VALUE(maxFractionNotConvergedCells, "[0")
        ATTRIBUTE_MAX_VALUE(maxFractionNotConvergedCells, "1]")
        ATTRIBUTE_DEFAULT_VALUE(maxFractionNotConvergedCells, "0.1")
        ATTRIBUTE_DISPLAYED_IF(maxFractionNotConvergedCells, "Level2")

        PROPERTY_DOUBLE(stabilityConvergenceThreshold, "convergence threshold on the converged-cell fraction")
        ATTRIBUTE_MIN_VALUE(stabilityConvergenceThreshold, "[0")
        ATTRIBUTE_MAX_VALUE(stabilityConvergenceThreshold, "1]")
        ATTRIBUTE_DEFAULT_VALUE(stabilityConvergenceThreshold, "0.005")
        ATTRIBUTE_DISPLAYED_IF(stabilityConvergenceThreshold, "Level2")

        PROPERTY_DOUBLE(maxChangeInGlobalIonizedH, "convergence threshold on the total ionized-H change")
        ATTRIBUTE_MIN_VALUE(maxChangeInGlobalIonizedH, "[0")
        ATTRIBUTE_MAX_VALUE(maxChangeInGlobalIonizedH, "1]")
        ATTRIBUTE_DEFAULT_VALUE(maxChangeInGlobalIonizedH, "0.001")
        ATTRIBUTE_DISPLAYED_IF(maxChangeInGlobalIonizedH, "Level2")

        // Diffuse reemission
        PROPERTY_DOUBLE(reemissionFraction, "fraction of the theoretical reemission probability that is realised")
        ATTRIBUTE_MIN_VALUE(reemissionFraction, "[0")
        ATTRIBUTE_MAX_VALUE(reemissionFraction, "1]")
        ATTRIBUTE_DEFAULT_VALUE(reemissionFraction, "1.0")
        ATTRIBUTE_DISPLAYED_IF(reemissionFraction, "Level3")

        PROPERTY_DOUBLE(transitionBlendWidth, "blending width in logU between the standard and transition tables")
        ATTRIBUTE_MIN_VALUE(transitionBlendWidth, "[0")
        ATTRIBUTE_MAX_VALUE(transitionBlendWidth, "1.5]")
        ATTRIBUTE_DEFAULT_VALUE(transitionBlendWidth, "0.3")
        ATTRIBUTE_DISPLAYED_IF(transitionBlendWidth, "Level2")

        PROPERTY_DOUBLE(maxHydrogenDensity, "maximum hydrogen number density ceiling in cm^-3 (0 means no ceiling)")
        ATTRIBUTE_MIN_VALUE(maxHydrogenDensity, "[0")
        ATTRIBUTE_DEFAULT_VALUE(maxHydrogenDensity, "100000")

        PROPERTY_DOUBLE(maxEmissionWavelength, "the upper limit of the gas emission wavelength grid")
        ATTRIBUTE_QUANTITY(maxEmissionWavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(maxEmissionWavelength, "1 micron")
        ATTRIBUTE_MAX_VALUE(maxEmissionWavelength, "1 m")
        ATTRIBUTE_DEFAULT_VALUE(maxEmissionWavelength, "1000 micron")
        ATTRIBUTE_DISPLAYED_IF(maxEmissionWavelength, "Level3")

        PROPERTY_INT(numEmissionWavelengths, "the number of bins in the gas emission wavelength grid")
        ATTRIBUTE_MIN_VALUE(numEmissionWavelengths, "50")
        ATTRIBUTE_MAX_VALUE(numEmissionWavelengths, "2000")
        ATTRIBUTE_DEFAULT_VALUE(numEmissionWavelengths, "250")
        ATTRIBUTE_DISPLAYED_IF(numEmissionWavelengths, "Level3")

        PROPERTY_BOOL(includeExtendedLines, "use Case B recombination tables and the atomic-model collisional "
                                            "solver")
        ATTRIBUTE_DEFAULT_VALUE(includeExtendedLines, "true")
        ATTRIBUTE_DISPLAYED_IF(includeExtendedLines, "Level3")

        PROPERTY_DOUBLE(lineLuminosityFloor, "the relative luminosity floor for extended-inventory lines")
        ATTRIBUTE_MIN_VALUE(lineLuminosityFloor, "[0")
        ATTRIBUTE_MAX_VALUE(lineLuminosityFloor, "1]")
        ATTRIBUTE_DEFAULT_VALUE(lineLuminosityFloor, "0")
        ATTRIBUTE_DISPLAYED_IF(lineLuminosityFloor, "Level3")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function opens the STAB tables for temperature and opacity (Standard and Transition),
        loads the reemission spectra table if diffuse reemission is enabled, initializes wavelength
        boundaries for the 5-bin calculation, builds the continuum emission wavelength grid,
        and initializes the PhotoIonizationSolver used for computing ion fractions at emission time. */
    void setupSelfBefore() override;

    //======== Capabilities =======

public:
    /** Returns true, indicating that this material mix has state variables beyond number density */
    bool hasExtraSpecificState() const override;

    /** Returns DynamicStateType::Primary, indicating that this material mix has a dynamic
        medium state with updates that are considered to affect primary emission. */
    DynamicStateType hasDynamicMediumState() const override;

    /** Returns true because emission is handled as continuous spectrum */
    bool hasContinuumEmission() const override;

    /** Returns true if this material mix supports scattering (for diffuse reemission). */
    bool hasScattering() const;

    /** This function returns true when diffuse reemission is turned on. */
    bool hasScatteringDispersion() const override;

    /** Returns a wavelength grid for the emission spectrum */
    DisjointWavelengthGrid* emissionWavelengthGrid() const override;

    //======== Medium state setup =======

public:
    /** Returns descriptors for snapshot-import columns required by this material mix.
        When abundanceMode=PerCell, requests 9 columns: 8 metal n_X/n_H ratios (C, N,
        O, Ne, Mg, Si, S, Fe) followed by n_He/n_H. SolarScaled requires no columns. */
    vector<SnapshotParameter> parameterInfo() const override;

    /** Returns descriptors for state variables used by this material mix. */
    vector<StateVariable> specificStateVariableInfo() const override;

    /** Initializes specific state variables from imported or default values. */
    void initializeSpecificState(MaterialState* state, double metallicity, double temperature,
                                 const Array& params) const override;

    //======== Medium state updates =======

    /** Updates ionization fractions, ionization parameter and 4-bin parameters based on the radiation field. */
    UpdateStatus updateSpecificState(MaterialState* state, const Array& Jv) const override;

    /** Determines if the medium state has converged based on ionization fractions and parameters. */
    bool isSpecificStateConverged(int numCells, int numUpdated, int numNotConverged, MaterialState* currentAggregate,
                                  MaterialState* previousAggregate) const override;

    //======== Low-level material properties =======

public:
    /** Returns the mean mass of a gas particle (weighted by composition). */
    double mass() const override;

    /** Returns a representative absorption cross section (in m2). */
    double sectionAbs(double lambda) const override;

    /** Returns the scattering cross section (zero for this material). */
    double sectionSca(double lambda) const override;

    /** Returns the total extinction cross section (equals absorption for this material). */
    double sectionExt(double lambda) const override;

    //======== High-level photon life cycle =======

    /** Returns the absorption opacity based on gas neutral fractions. */
    double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** Returns scattering opacity for diffuse reemission if enabled. */
    double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** Returns the extinction opacity (equals absorption for this material). */
    double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** Performs peel-off scattering for diffuse reemission. Returns false: no
        fluorescence-as-secondary-emission. */
    bool peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs, Direction bfky,
                           const MaterialState* state, const PhotonPacket* pp) const override;

    /** Performs scattering event including possible diffuse reemission. */
    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======== Secondary emission =======

    /** Returns the continuum emission spectrum (free-bound, free-free, two-photon). */
    Array emissionSpectrum(const MaterialState* state, const Array& Jv) const override;

    /** Returns true: DIG emits discrete emission lines computed from PIO-style ion balance. */
    bool hasLineEmission() const override;

    /** Returns the rest-frame wavelength centers of all emission lines. */
    Array lineEmissionCenters() const override;

    /** Returns the particle masses for thermal line broadening. */
    Array lineEmissionMasses() const override;

    /** Returns the per-line luminosities [W] computed from T, ne, and ion fractions. */
    Array lineEmissionSpectrum(const MaterialState* state, const Array& Jv) const override;

    /** Returns the fixed temperature used for this material. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //======== 5-Bin Radiation Field Characterization =======

private:
    // update temperature from STAB table if enabled
    void updateTemperatureFromStab(MaterialState* state) const;

    // Calculate 5-bin ratio parameters R2, R3, R4, and R5 from radiation field
    void calculate5BinRatioParameters(const Array& Jv, double& logR2, double& logR3, double& logR4,
                                      double& logR5) const;

    // Calculate average intensity in energy bins
    void calculateBinAverages(const Array& Jv, double* binAverages) const;

    // Convert Rydberg energy to wavelength in metres
    double rydbergToWavelength(double energy_ryd) const;

    // Per-channel reemission probabilities for a cell/wavelength; defined in the .cpp along with the
    // channel enumeration and the probabilities array size, none of which is needed outside this
    // class's own translation unit
    struct ReemissionData;

    // Reemission calculation methods (adapted from CMACIONIZE)
    void calculateReemissionProbabilities(const MaterialState* state, double lambda, ReemissionData& data) const;
    int selectReemissionChannel(const MaterialState* state, double lambda, const PhotonPacket* pp) const;
    double reemissionChannelWeights(const MaterialState* state, const ReemissionData& data, double weights[]) const;
    double sampleFromTemperatureDependentSpectrum(int spectrumType, double temperature, Random* random) const;
    void interpolateReemissionSpectrum(int spectrumType, double temperature, vector<double>& wavelengths,
                                       vector<double>& cumulativeDist) const;
    double sampleHydrogenLymanContinuum(double temperature) const;
    double sampleHeliumLymanContinuum(double temperature) const;
    double sampleHeliumTwoPhotonContinuum(double temperature) const;

    // Hydrogen and Helium photoionization cross-section using Verner+ 96
    double getHydrogenCrossSection(double frequency) const;
    double getHeliumCrossSection(double frequency) const;

    // Each cell stores the extinction opacity at the ionizing grid points (all other grid points hold zero opacity by
    // construction) and the effective re-emission fraction at the grid points of the re-emission range, as
    // consecutive custom state variables; the absorption and scattering opacities follow from these two
    int _indexFirstOpacityExt = 0;          // first index of the extinction opacity array
    int _indexFirstReemissionFraction = 0;  // first index of the effective re-emission fraction array
    int _numIonizingOpacityPoints = 0;      // grid points 0 .. n-1 have lambda <= _lambdaLow
    int _firstReemissionOpacityPoint = 0;   // first grid point with _lambdaBin5 <= lambda <= _lambdaH
    int _numReemissionOpacityPoints = 0;    // zero when reemissionFraction() is zero

    // Radiation field wavelength grid caching
    Array _opacityWavelengthGrid;       // Cached wavelength grid for opacity calculations
    double _opacityMinWavelength = 0.;  // left border of the first radiation field bin

    // Pre-compute opacity arrays for the radiation field wavelength grid
    void precomputeOpacityArrays(MaterialState* state, const double* solverIonFracs) const;

    // Get the interpolated opacity from pre-computed state variables
    // opacityType: 0=absorption, 1=scattering, 2=extinction
    double interpolateOpacityFromState(double lambda, MaterialState* state, int opacityType) const;

    // Initialize the opacity wavelength grid from the radiation field (called during setup)
    void initializeOpacityWavelengthGrid();

    // Radiation field characterization
    double calculateIonizationParameter(const Array& Jv, double nH) const;

    // convert total to Hydrogen number density
    double convertTotalDensityToHydrogenDensity(double n_total, double He_abundance) const;

    // Build the per-cell 8-element metal abundance vector {n_C, n_N, n_O, n_Ne, n_Mg, n_Si, n_S, n_Fe}
    // (number ratios n_X/n_H) for the analytical PIO/emission solver. In SolarScaled mode the vector
    // is computed analytically from solar pattern * Zfrac * g * scalings; in PerCell mode it is
    // read from custom state slots that were filled at init from snapshot params.
    void buildPerCellAbundances(const MaterialState* state, double abund[8]) const;

    // Resolve a cell's (dN, dC) abundance perturbation relative to solar*Zfrac. SolarScaled mode
    // returns log10 of solarScalingN/C; PerCell mode derives from the imported (n_N, n_C / n_H)
    // and the cell metallicity.
    void computeCellDeltas(const MaterialState* state, double& dN, double& dC) const;

    // Bracket a query value on a sorted delta-axis (axis_values, axis_deltaIds). Endpoint-clamps
    // outside the cross. Returns the two delta_id indices into the stab and the linear weight w
    // such that  reconstructed(value) = (1-w)*tab(deltaIdLo) + w*tab(deltaIdHi).
    void bracketDeltaAxis(double value, const vector<double>& axis_values, const vector<int>& axis_deltaIds,
                          int& deltaIdLo, int& deltaIdHi, double& w) const;

    // Dual-table system helpers
    enum class TableSelection { Standard, Transition, Blend };
    TableSelection selectTable(double logU) const;  // Selects which table to use based on logU
    double getBlendingWeight(double logU) const;    // Returns weight for transition table (0-1)

    // Data members - Dual tables for Standard and Transition regions (multiZ axes).
    StoredTable<8> _standardTemperatureTable;    // 8D Standard temperature: Z,n_H,logU,logR2..5,delta_id
    StoredTable<7> _transitionTemperatureTable;  // 7D Transition temperature: Z,n_H,logU,logR2..5
    StoredTable<9> _standardOpacityTable;        // 9D Standard opacity: lambda,logU,logR2..5,Z,n_H,delta_id
    StoredTable<8> _transitionOpacityTable;      // 8D Transition opacity: lambda,logU,logR2..5,Z,n_H
    StoredTable<3> _reemissionSpectraTable;      // 3D cumulative probability table

    // Delta-map table for the multiZ stab (loaded whenever any temperature/opacity stab is on).
    // One StoredTable instance per quantity; both share the same underlying file resource.
    StoredTable<1> _deltaNMapTable;  // 1D map: delta_id -> delta_N (dex)
    StoredTable<1> _deltaCMapTable;  // 1D map: delta_id -> delta_C (dex)
    // Stab dN/dC cross axes for separable bracket+linterp reconstruction (built at setup).
    // _dNAxisValues is sorted; _dNAxisDeltaIds holds the corresponding delta_id index in the stab.
    vector<double> _dNAxisValues;
    vector<int> _dNAxisDeltaIds;
    vector<double> _dCAxisValues;
    vector<int> _dCAxisDeltaIds;
    int _deltaIdCentre = -1;  // delta_id of the (dN=0, dC=0) centre slice
    double _stabDNmin = 0.;
    double _stabDNmax = 0.;
    double _stabDCmin = 0.;
    double _stabDCmax = 0.;

    // Wavelength boundaries for different ionization thresholds
    double _lambdaH = 0;   // Hydrogen ionization threshold wavelength (1 Ryd)
    double _lambdaHe = 0;  // Helium ionization threshold wavelength (~1.8 Ryd)

    // 5-bin wavelength boundaries (calculated from energy boundaries)
    double _lambdaBin1 = 0;  // 1.00 Ryd boundary
    double _lambdaBin2 = 0;  // 1.80 Ryd boundary
    double _lambdaBin3 = 0;  // 2.58 Ryd boundary
    double _lambdaBin4 = 0;  // 3.52 Ryd boundary
    double _lambdaBin5 = 0;  // 4.00 Ryd boundary
    double _lambdaBin6 = 0;  // 6.00 Ryd boundary

    // Radiation field bins assigned to each of the five energy bins (by characteristic wavelength), and the summed
    // width of those bins, set up once from the radiation field wavelength grid
    std::array<vector<int>, 5> _energyBinRFBins;
    std::array<double, 5> _energyBinRFWidth{};

    // Wavelength bounds
    double _lambdaLow = 0;   // Low energy cutoff @ 912 A
    double _lambdaHigh = 0;  // High energy cutoff @ 152 A

    // Emission wavelength grid
    DisjointWavelengthGrid* _emissionWavelengthGrid = nullptr;

    // Cached density ceiling in m^-3 (0 = disabled), set in setupSelfBefore from maxHydrogenDensity
    double _maxNH_m3 = 0.;

    // Cached upper edge of the transition table's logU axis, where it meets the standard
    // table, set in setupSelfBefore. Used when selecting between the two tables.
    double _transitionLogUMax = 0.;

    // Custom state variable indices
    int _indexHeliumAbundance = 0;
    int _indexHNeutralFraction = 0;
    int _indexHeNeutralFraction = 0;
    int _indexIonizationParameter = 0;
    int _indexLogR2 = 0;                 // log10(R2) ratio parameter
    int _indexLogR3 = 0;                 // log10(R3) ratio parameter
    int _indexLogR4 = 0;                 // log10(R4) ratio parameter
    int _indexLogR5 = 0;                 // log10(R5) ratio parameter
    int _indexNHIonized = 0;             // ionized hydrogen number density (for global convergence criterion)
    int _indexFirstIonFractionDiag = 0;  // starting index for diagnostic ion fraction block (per-cell, dimensionless)
    int _indexFirstIonFractionAgg =
        0;  // starting index for n_HII-weighted ion fraction aggregates (numbervolumedensity)
    int _indexFirstMetalAbund =
        -1;  // starting index for 8 per-cell metal n_X/n_H slots (only when abundanceMode=PerCell; -1 otherwise)

    // Cached lower density bound for opacity scaling (set in setupSelfBefore)
    double _nHMinOpacity = 0;

    // PIO-style emission: solver for computing ion fractions from Jv + T
    PhotoIonizationSolver _emissionSolver;

    // Line inventory, Case B tables and atomic models shared by the line-emission calculations
    // below; owning this instance (rather than sharing process-global state) means its resources
    // are released when this DiffuseIonizedGasMix instance is, not at unspecified static
    // destruction order
    GasLineEmission _gasLineEmission;

    // Line emission data (initialized in setupSelfBefore)
    int _numLines = 0;
    Array _lineCenters;
    Array _lineMasses;
    // registry indices of the lines this mix emits: all built-in lines plus the
    // extended-inventory lines that fall inside the emission wavelength grid
    vector<int> _activeLines;
    // solve groups: active-line indices per atomic-model slot (one level-population solve per
    // group per cell), with each line's transition index; lines without a model use the
    // per-line path
    vector<int> _groupSlots;
    vector<vector<int>> _groupLines;
    vector<int> _lineTransitions;
    vector<int> _lineGroup;  // group index per active line, -1 = per-line path

    // Convergence stability tracking
    mutable vector<double> _convergedFractionHistory;
    mutable size_t _convergenceHistorySize = 3;
};

//////////////////////////////////////////////////////////////////////

#endif  // DIFFUSEIONIZEDGASMIX_HPP
