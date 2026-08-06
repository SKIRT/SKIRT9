/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PHOTOIONIZATIONSOLVER_HPP
#define PHOTOIONIZATIONSOLVER_HPP

#include "Array.hpp"

//////////////////////////////////////////////////////////////////////

/** The PhotoIonizationSolver class implements a first-principles photoionization solver for
    computing the ionization state, equilibrium temperature, opacity, and line emissivities of
    ionized gas given a radiation field J_lambda (mean intensity per unit wavelength).
    It handles H, He, C, N, O, Ne, Mg, Si, S, and Fe through a unified pipeline.

    Note: in SKIRT, the array conventionally called "Jv" stores J_lambda [W m^-3 sr^-1],
    not J_nu. The solver converts from wavelength to frequency/energy internally.

    The solver implements:
    - Ionization balance via the Katz, Weinberg & Hernquist (1996) recursion
    - Temperature from explicit heating = cooling equilibrium (bisection)
    - Heating from full J_lambda integration (no 5-bin compression)
    - Cooling from analytic H/He channels plus pre-tabulated metal line cooling
    - Opacity computed analytically from ion fractions and Verner+ cross-sections

    The solver is designed as a standalone utility class (not a SimulationItem) so it can be
    embedded as an inline ion-balance helper by any material mix that needs first-principles
    photoionisation. The interface accepts a prior state (T_prior, ion fractions) as the
    initial estimate, which accelerates convergence. All internal computations use CGS units
    (cm^-3 for densities, erg/s/cm^3 for rates, eV for energies). */
class PhotoIonizationSolver final
{
public:
    //======== Element indexing =======

    /** Number of elements tracked by the solver. */
    static constexpr int numElements = 10;

    /** Number of ionization stages per element (includes fully ionized stage). */
    static constexpr int numStages[] = {2, 3, 7, 8, 9, 9, 11, 13, 5, 17};
    //                                  H  He  C   N   O  Ne  Mg  Si  S  Fe

    /** Total number of ion fractions stored (sum of numStages). */
    static constexpr int totalStages = 2 + 3 + 7 + 8 + 9 + 9 + 11 + 13 + 5 + 17;  // = 84

    /** Index of first ionization stage for each element in the ion fraction array. */
    static constexpr int stageOffset[] = {0, 2, 5, 12, 20, 29, 38, 49, 62, 67};

    //======== Result structure =======

    /** Result of solving the ionization and thermal equilibrium for a single cell. */
    struct CellResult
    {
        double Teq = 1e4;              ///< equilibrium temperature [K]
        double ne = 0.;                ///< electron number density [cm^-3]
        double heatingRate = 0.;       ///< photoheating rate [erg s^-1 cm^-3]
        double coolingRate = 0.;       ///< cooling rate [erg s^-1 cm^-3]
        double ionFracs[totalStages];  ///< ion fractions for all species (sum per element = 1)
    };

    //======== Construction =======

    /** Default constructor. */
    PhotoIonizationSolver() = default;

    /** Initializes the solver for a given radiation field wavelength grid.
        \param lambdav wavelength grid bin centers [m], must be sorted ascending
        \param dlambdav wavelength grid bin widths [m]
        \param coolingTable true to load cooling table, false to skip */
    void initialize(const Array& lambdav, const Array& dlambdav, bool coolingTable = false);

    //======== Main solve interface =======

    /** Solves the ionization balance and thermal equilibrium for a single cell.
        \param Jv mean intensity J_lambda at each wavelength bin [W m^-3 sr^-1]
        \param nH hydrogen number density [cm^-3]
        \param yHe helium abundance by number relative to H
        \param metalAbundances number abundances relative to H for {C, N, O, Ne, Mg, Si, S, Fe}
        \param Tprior prior temperature estimate [K]
        \param ionFracsPrior prior ion fractions (can be nullptr for default initialization)
        \param TdampFactor maximum factor by which T can change per iteration (0 or <=1 = no damping) */
    CellResult solve(const Array& Jv, double nH, double yHe, const double metalAbundances[8], double Tprior,
                     const double* ionFracsPrior, double TdampFactor = 0.) const;

    /** Solves the ionization balance at a fixed temperature. The temperature is kept
        fixed; only the ionization balance is solved (no T bisection). Photoionization rates
        from Jv are included, so this is a photo+collisional hybrid. Heating and cooling rates
        are computed for diagnostics but do not affect T. Used for CIE cells (shock-heated gas
        retaining snapshot T) and for re-solving ion fractions after T sign-flip damping.
        \param Jv mean intensity J_lambda at each wavelength bin [W m^-3 sr^-1]
        \param nH hydrogen number density [cm^-3]
        \param yHe helium abundance by number relative to H
        \param metalAbundances number abundances relative to H for {C, N, O, Ne, Mg, Si, S, Fe}
        \param T fixed temperature [K]
        \param ionFracsPrior prior ion fractions (can be nullptr for default initialization)
        \param useEffectiveStages if true, use maxIonizationEnergy-capped stages (for sign-flip
               re-solve consistency); if false, use full stages (for CIE cells) */
    CellResult solveIonizationAtFixedT(const Array& Jv, double nH, double yHe, const double metalAbundances[8],
                                       double T, const double* ionFracsPrior, bool useEffectiveStages = false) const;

    //======== Opacity interface =======

    /** Returns the absorption opacity [cm^-1] at wavelength lambda [m] given the ion fractions
        and hydrogen number density. */
    double opacityAbs(double lambda, const double ionFracs[totalStages], double nH, double yHe,
                      const double metalAbundances[8]) const;

    /** Initializes ion fractions to a default state (fully ionized for H/He, singly ionized
        for metals). */
    void initializeIonFractions(double ionFracs[]) const;

    /** Sets the maximum ionization energy threshold [eV]. Ion stages with ionization potential
        above this threshold are excluded from the PIO solve (capped). CIE mode always uses
        full stages. Negative value means no cap (all stages active). Must be called after
        initialize(). */
    void setMaxIonizationEnergy(double maxEV);

    /** Enables or disables the cosmic-ray background heating term (Indriolo+ 2007; Wolfire+ 1995).
        When enabled, adds a heating rate of scaleFactor * zeta_CR * 35 eV * nHI per cell.
        The default CR ionization rate is zeta_CR = 2e-16 s^-1 (Indriolo+ 2007). */
    void setCosmicRayHeating(bool enabled, double scaleFactor = 1.0);

    /** Returns the total number of effective (active) ion stages across all elements. */
    int totalEffectiveStages() const;

    /** Returns the effective number of stages for element e. */
    int effectiveStages(int e) const { return _nStagesEff[e]; }

    /** Returns the number of wavelength bins in the solver's grid. */
    int numBins() const { return _numBins; }

    /** Computes photoionization rates [s^-1] for all ions by integrating 4*pi*J_nu*sigma(nu)/(h*nu)
        over the wavelength grid. Results stored in gamma[]. */
    void computePhotoionizationRates(const Array& Jv, double gamma[]) const;

private:
    //======== Internal methods =======

    /** Computes the photoheating rate [erg s^-1 cm^-3] from the radiation field and ion fractions.
        Integrates 4*pi*J_nu*sigma(nu)*(h*nu - IP)/(h*nu) over the wavelength grid, weighted by
        ion number densities n_ion. The effNS[] array specifies effective stages per element. */
    double computeHeatingRate(const Array& Jv, const double ionFracs[], double nH, double yHe,
                              const double metalAbundances[8], double ne, const int effNS[]) const;

    /** Computes the total cooling rate [erg s^-1 cm^-3] from analytic H/He cooling channels
        (recombination, free-free, collisional ionization, collisional excitation) plus metal
        line cooling. The effNS[] array specifies effective stages per element. */
    double computeCoolingRate(double T, const double ionFracs[], double nH, double yHe, const double metalAbundances[8],
                              double ne, const int effNS[]) const;

    /** Solves the ionization balance for all species at temperature T and electron density ne,
        given photoionization rates gamma[]. Uses the Katz+ (1996) recursion for metals and
        iterative solution for H/He. Updates ionFracs[] and returns the new electron density.
        The effNS[] array specifies effective stages per element. */
    double solveIonizationBalance(double T, double ne, const double gamma[], double nH, double yHe,
                                  const double metalAbundances[8], double ionFracs[], const int effNS[]) const;

    /** Loads the CHIANTI metal line cooling table from file.
        Parses ion names to build the mapping from ionFracs index to table quantity index. */
    void loadCoolingTable();

    /** Returns the metal line cooling rate [erg/s per ion] for the given ionFracs index
        by bilinear interpolation in log T and log ne from the pre-loaded table.
        Returns 0 if no table is loaded or the ion has no cooling data. */
    double interpolateMetalCooling(int ionFracsIdx, double T, double ne) const;

    //======== Wavelength grid data =======

    int _numBins = 0;                ///< number of wavelength bins
    std::vector<double> _lambda;     ///< wavelength bin centers [m]
    std::vector<double> _dlambda;    ///< wavelength bin widths [m]
    std::vector<double> _energy;     ///< photon energy per bin [eV]
    std::vector<double> _energyErg;  ///< photon energy per bin [erg]

    //======== Effective stage counts (may be capped by maxIonizationEnergy) =======

    int _nStagesEff[numElements] = {2, 3,  7,  8, 9,
                                    9, 11, 13, 5, 17};  ///< effective stages per element (default = full)

    //======== Metal cooling table data =======

    static constexpr int _metalOffset = 5;                            ///< first metal ionFracs index (CI)
    static constexpr int _numMetalIons = totalStages - _metalOffset;  ///< 55 metal ion stages

    bool _hasCoolingTable = false;    ///< whether a cooling table is loaded
    int _coolNT = 0;                  ///< number of T grid points
    int _coolNNe = 0;                 ///< number of ne grid points
    std::vector<double> _coolLogT;    ///< log10(T/K) grid
    std::vector<double> _coolLogNe;   ///< log10(ne/cm^-3) grid
    std::vector<double> _coolData;    ///< cooling[ion][iT][iNe] flattened, [erg/s per ion]
    std::vector<int> _stabIdxForIon;  ///< maps ionFracs index -> .stab quantity index (-1 if none)

    //======== Cosmic-ray background heating =======

    bool _cosmicRayHeating = false;  ///< whether to include CR background heating
    double _cosmicRayScale = 1.0;    ///< multiplicative scale factor for CR heating rate
};

//////////////////////////////////////////////////////////////////////

#endif
