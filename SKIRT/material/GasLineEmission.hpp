/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GASLINEEMISSION_HPP
#define GASLINEEMISSION_HPP

#include "StoredTable.hpp"
class Log;
class SimulationItem;

//////////////////////////////////////////////////////////////////////

/** GasLineEmission is a helper class that provides species-agnostic building blocks for gas
    line emission, offered for use by material mixes:

    - Line inventory: 26 built-in lines (H and He recombination lines, optical forbidden metal
      lines) indexed by LineIndex, plus the extended inventory appended at setup by
      initializeExtendedLineRegistry().
    - Recombination lines (Case B): L = eps(T, n_e) n_e n_ion V from emissivity tables (Storey &
      Sochi 2015 for H I, Porter et al. 2012 for He I, Storey & Hummer 1995 for He II) loaded by
      initializeRecombinationTables(); without tables the H lines use the legacy P_B form (Storey
      & Hummer 1995; Hui & Gnedin 1997; McClymont, Smith & Tacchella 2025) and the He lines are
      zero. Lyman-alpha always uses the legacy form (the Case B tables exclude the Lyman series).
    - Collisional lines in the nebular limit (level populations set by electron collisions at the
      local T and n_e, no radiative pumping): the statistical-equilibrium solver on atomic models
      loaded by initializeAtomicModels(), or the legacy precomputed q_col(T, n_e) tables if this
      function has not been called.
    - Statistical equilibrium: the general level-population solver (collisions with any set of
      partners, spontaneous decay, optional radiative pumping) on atomic models read from the
      NonLTELineGasMix species files. The solver follows Kosei Matsumoto's NonLTELineGasMix
      implementation (Matsumoto et al. 2023) and is shared by NonLTELineGasMix and
      DiffuseIonizedGasMix.

    Clients should construct a GasLineEmission instance as data member and call initialize()
    before calling any other functions on the instance. All initilization must be performed
    during setup in the master thread. Once the instance is initialized, all other functions\
    are thread-safe.
*/
class GasLineEmission final
{
    // ============== Constructing and destructing ==============

public:
    /** The constructor initializes the instance to an invalid state. Call initialize() before
        calling any other functions. Failing to do so causes undefined behavior. */
    GasLineEmission();

    /** The destructor releases any resources acquired by this class instance. */
    ~GasLineEmission();

    /** The copy constructor and copy assignment operator are deleted because an instance may own
        resources (such as stored tables) that cannot be meaningfully shared by naive member-wise
        copying. */
    GasLineEmission(const GasLineEmission&) = delete;
    GasLineEmission& operator=(const GasLineEmission&) = delete;

    // ============== Public constants and data types ==============

public:
    /** Line indices for the emission line array. */
    enum LineIndex : int {
        // H recombination lines (Case B)
        Lya = 0,   // Lyman-alpha   1215.67 A
        Ha,        // Balmer-alpha   6562.80 A
        Hb,        // Balmer-beta    4861.33 A
        Hg,        // Balmer-gamma   4340.46 A
        Hd,        // Balmer-delta   4101.73 A
        HeBalmer,  // Balmer-epsilon 3970.07 A
        Paa,       // Paschen-alpha  18751 A
        Pab,       // Paschen-beta   12818 A
        Bra,       // Brackett-alpha 40512 A

        // He recombination lines (Case B)
        HeI5876,   // He I  5876 A
        HeI6678,   // He I  6678 A
        HeI7065,   // He I  7065 A
        HeI10830,  // He I  10830 A
        HeII1640,  // He II 1640 A
        HeII4686,  // He II 4686 A

        // Collisional lines (optical forbidden metal lines)
        NII6548,   // [NII] 6548
        NII6583,   // [NII] 6583
        OI6300,    // [OI] 6300
        OI6364,    // [OI] 6364
        OII3729,   // [OII] 3729  (2D 5/2)
        OII3726,   // [OII] 3726  (2D 3/2)
        OIII4363,  // [OIII] 4363
        OIII4959,  // [OIII] 4959
        OIII5007,  // [OIII] 5007
        SII6716,   // [SII] 6716
        SII6731,   // [SII] 6731
        numLines
    };

    /** Record representing one registry line: carrierIonIndex/elementIndex address the consumer's
        ion-fraction and abundance arrays (PhotoIonizationSolver layout); both are -1 for
        recombination lines, whose recombining ion (H II, He II, He III) is selected by the family
        flags. */
    struct LineDef
    {
        double wavelength;    // rest-frame wavelength [m]
        double mass;          // particle mass for Doppler broadening [kg]
        int carrierIonIndex;  // PhotoIonizationSolver ion stage index, -1 for recombination lines
        int elementIndex;     // abundance element index (0=C..7=Fe), -1 for recombination lines
        bool isHeIRecomb;     // He I recombination line (recombining ion He II)
        bool isHeIIRecomb;    // He II recombination line (recombining ion He III)
    };

    /** Record representing one atomic species for the extended inventory: resource base name (e.g.
        "N_II") plus the consumer's carrier indices. */
    struct SpeciesSpec
    {
        string name;          // atomic data file base name
        int carrierIonIndex;  // ion stage index in the consumer's ion fraction array
        int elementIndex;     // abundance element index
    };

    /** Record holding collisional transition data for one collision partner. */
    struct ColPartner
    {
        string name;                 // partner name ("e-", "H", ...)
        vector<double> T;            // temperature grid [K]
        vector<int> indexUpCol;      // upper level index per collisional transition
        vector<int> indexLowCol;     // lower level index per collisional transition
        vector<vector<double>> Kul;  // de-excitation rate [m3/s] on the T grid
    };

    /** Record holding atomic data for one species in SI units, as read by loadAtomicModel(). */
    struct AtomicModel
    {
        double mass{0.};                // particle mass [kg]
        vector<double> energy;          // level energies [J]
        vector<double> weight;          // statistical weights
        vector<int> indexUpRad;         // upper level index per radiative transition
        vector<int> indexLowRad;        // lower level index per radiative transition
        vector<double> einsteinA;       // Einstein A coefficient [1/s]
        vector<double> einsteinBul;     // Einstein B_ul (per-wavelength convention)
        vector<double> einsteinBlu;     // Einstein B_lu (per-wavelength convention)
        vector<double> center;          // line center wavelength [m]
        vector<double> branchRatio;     // A divided by the sum of A from the same upper level
        vector<ColPartner> colPartner;  // collision partners

        int numLevels() const { return energy.size(); }
        int numLines() const { return einsteinA.size(); }
        int numColPartners() const { return colPartner.size(); }
    };

    /** Per-cell inputs for solveLevelPopulations(). An empty meanJ skips radiative pumping. */
    struct Environment
    {
        double Tkin{0.};          // kinetic temperature [K]
        double nTotal{0.};        // species number density [m^-3]
        vector<double> nPartner;  // collision partner densities [m^-3]
        vector<double> meanJ;     // per-line mean intensity; empty -> pumping skipped
    };

    // ============== Initialization -- call in main thread ==============

public:
    /** Initializes the line registry with the built-in lines (the first numLines entries, in
        LineIndex order). The \em item argument is remembered for use throughout the class to
        retrieve the logger and the units system. This function must be called before calling any
        other function. */
    void initialize(const SimulationItem* item);

    /** Loads the Case B emissivity tables from resources as memory mapped stored tables. Does
        nothing after being called once. */
    void initializeRecombinationTables();

    /** Loads the atomic models of the built-in collisional lines' carrier species and maps each
        line to its transition by nearest wavelength. Does nothing after being called once. */
    void initializeAtomicModels();

    /** Calls initializeRecombinationTables() and initializeAtomicModels() and then appends every
        recombination table listed in the specified per-set wavelength index files and every
        radiative transition of each species loadable from resources. Absent species are silently
        skipped and built-in lines are not duplicated. Does nothing after being called once. */
    void initializeExtendedLineRegistry(const vector<SpeciesSpec>& species);

    // ============== Recombination lines -- thread-safe after initialization ==============

public:
    /** Returns a reference to the line registry: the first numLines entries are the built-in lines
        (LineIndex order), followed by any lines added by initializeExtendedLineRegistry(). */
    const vector<LineDef>& lineRegistry() const;

    /** Returns legacy H line luminosity [W] (Lya through Bra): h nu P_B(T, n_e) gammaHI nHI V,
        with P_B the Case B probability that a recombination emits the line. Densities in cm^-3,
        volume in cm^3. */
    double hydrogenLineLuminosity(int lineIdx, double T, double ne, double gammaHI, double nHI, double V_cm3) const;

    /** Returns H or He recombination line luminosity [W] (Lya through HeII4686): eps(T, ne) ne
        nIon V from the emissivity tables when loaded, else the legacy form for H (which uses
        gammaHI, nHI) and zero for He. nIon is the recombining ion density (H II, He II or He III).
        Densities in cm^-3, volume in cm^3. */
    double recombinationLineLuminosity(int lineIdx, double T, double ne, double nIon, double gammaHI, double nHI,
                                       double V_cm3) const;

    // ============== Collisionally excited lines -- thread-safe after initialization ==============

public:
    /** Returns nebular-limit line luminosity [W] for a collisional line (NII6548 through SII6731,
        or an extended-inventory line): level populations from electron collisions at (T, ne)
        without radiative pumping, times nIon V. Uses the atomic models when loaded, else the
        legacy q_col tables. Densities in cm^-3, volume in cm^3. */
    double collisionalLineLuminosity(int lineIdx, double T, double ne, double nIon, double V_cm3) const;

    /** Returns the model slot of the species carrying the given registry line, or -1 if the line
        is not served by a loaded atomic model. Lines with the same slot share one solve. */
    int lineModelSlot(int lineIdx) const;

    /** Returns the transition index of the given registry line within its model (valid when
        lineModelSlot() >= 0). */
    int lineTransition(int lineIdx) const;

    // ============== Statistical equilibrium -- thread-safe after initialization ==============

public:
    /** Fills the atomic model from the resource files NAME_Mass.txt, NAME_Energy.txt,
        NAME_Rad_Coeff.txt and NAME_Col_PARTNER_Temp.txt / NAME_Col_PARTNER_Coeff.txt (NAME the
        species name, PARTNER each collision partner), keeping at most maxNumLevels levels. The
        resource files are opened silently, so a caller should issue its own summary log message
        for the operation as a whole after loading a species or a group of species. */
    void loadAtomicModel(const string& speciesName, const vector<string>& partnerNames, int maxNumLevels,
                         AtomicModel& model) const;

    /** Returns the atomic model loaded by initializeAtomicModels() into the given slot (see
        lineModelSlot()). */
    const AtomicModel& atomicModel(int slot) const;

    /** Solves the statistical-equilibrium rate matrix and returns level populations [m^-3]
        normalized to env.nTotal. Throws FatalError if the matrix is singular or the
        solution is not finite. */
    vector<double> solveLevelPopulations(const AtomicModel& model, const Environment& env) const;

    /** Returns the line power densities [W m^-3] for all radiative transitions: n_up A h c /
        lambda. */
    vector<double> lineEmissivities(const AtomicModel& model, const vector<double>& pops) const;

    // ================== Private data types and data members ==================

private:
    // the simulation item passed to initialize() and the logger retrieved from it
    const SimulationItem* _item{nullptr};
    Log* _log{nullptr};

    // the mutable line registry:
    //  - initialize() sets the built-in lines
    //  - initializeExtendedLineRegistry() appends to it in place
    vector<LineDef> _registry;

    // registry mapping recombination lines to loaded Case B emissivity tables; filled by
    // initializeRecombinationTables() and extended by initializeExtendedLineRegistry().
    struct RecombRegistry
    {
        bool ready = false;

        // one aggregated Case B emissivity cube Emis(line, T, ne) per species: [0] H I, [1] He I, [2] He II
        StoredTable<3> cube[3];

        // transition metadata read from each species' line index map (index is 1-based, matching the cube line axis)
        struct MapRow
        {
            int index;
            int upper;
            int lower;
            double wav_m;
        };
        vector<MapRow> map[3];

        // per registry line: which cube serves it and at which 1-based line index (cubeId -1 = not a Case B line)
        struct Ref
        {
            int cubeId = -1;
            int lineIdx = 0;
        };
        vector<Ref> table;

        bool loaded(int idx) const { return table[idx].cubeId >= 0; }
    };
    RecombRegistry _recombRegistry;

    // maps the built-in collisional lines (and any added by initializeExtendedLineRegistry())
    // to loaded atomic models and transitions; filled by initializeAtomicModels() and extended by
    // initializeExtendedLineRegistry()
    struct AtomicLineRegistry
    {
        bool ready = false;
        vector<AtomicModel> models;  // one per carrier species
        vector<string> modelNames;   // species name per model slot
        vector<int> lineModel;       // model slot per line, -1 = none
        vector<int> lineTransition;  // transition index within the model
    };
    AtomicLineRegistry _atomicRegistry;
};

//////////////////////////////////////////////////////////////////////

#endif
