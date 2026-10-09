/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DiffuseIonizedGasMix.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "GasContinuumEmission.hpp"
#include "GasLineEmission.hpp"
#include "Log.hpp"
#include "MaterialState.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"
#include "SnapshotParameter.hpp"
#include "StringUtils.hpp"
#include "VernerCrossSections.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // Minimum ionization parameter for emission: cells with logU <= this get T = _minTEmit
    // in updateTemperatureFromStab() and are also directly gated in emissionSpectrum().
    static constexpr double _minLogUEmit = -6.50;

    // Temperature floor [K] applied to cells below the logU emission gate
    static constexpr double _minTEmit = 300.0;

    // Log-space floor for logU and logR2..logR5: uninitialised cells sit here, and computed
    // values are clamped here. "At the floor" is detected with `value <= _logFloor`.
    static constexpr double _logFloor = -99.0;

    // Solar reference metallicity (GASS10).
    static constexpr double _solarZ = 0.014;

    // Solar abundances (GASS10) for Z-scaling of metal emission.
    // Order: C, N, O, Ne, Mg, Si, S, Fe.
    static constexpr double _solarAbundances[8] = {2.69e-4, 6.76e-5, 4.90e-4, 8.51e-5,
                                                   3.98e-5, 3.24e-5, 1.32e-5, 3.16e-5};

    // Per-element gas-phase depletion factors following Gunasekera et al. (2023): the
    // Jenkins (2009) depletion pattern (Table 4) evaluated at F* = 0.5, with no depletion
    // applied to sulphur. These multiply the solar abundance pattern in
    // buildPerCellAbundances so that the analytical ionization and line-emission step
    // uses the same gas-phase composition that Cloudy used when the temperature and
    // opacity tables were generated. Order matches _solarAbundances.
    static constexpr double _gunasekeraDepletion[8] = {
        0.688023,  // C  : J09 F*=0.5
        0.778037,  // N  : J09 F*=0.5
        0.753442,  // O  : J09 F*=0.5
        1.000000,  // Ne : noble gas, no depletion
        0.170179,  // Mg : J09 F*=0.5
        0.161614,  // Si : J09 F*=0.5
        1.000000,  // S  : Gunasekera (2023) leaves S undepleted
        0.025471,  // Fe : J09 F*=0.5
    };

    // Constants from CMACIONIZE (same values as in PhysicalDiffuseReemissionHandler)
    static constexpr double helium19p8EvFreq = 4.788e15;
    static constexpr double heliumTpcHIonizingFraction = 0.56;
    static constexpr double heliumLyaOtsCoeff = 77.0;

    // 4-bin energy structure boundaries (in Rydberg)
    static constexpr double energyBin1Ryd = 1.0;
    static constexpr double energyBin2Ryd = 1.8;
    static constexpr double energyBin3Ryd = 2.58;
    static constexpr double energyBin4Ryd = 3.52;
    static constexpr double energyBin5Ryd = 4.0;
    static constexpr double energyBin6Ryd = 6.0;

    // Define Rydberg energy (1 Rydberg = 13.6057 eV)
    constexpr double rydbergEnergyEv = 13.6057;

    // Numerical safety thresholds
    constexpr double minDivisionSafe = 1.e-12;

    // Physical bounds for neutral fractions
    constexpr double minNeutralFraction = 1.e-8;
    constexpr double maxNeutralFraction = 1.0;

    // Diagnostic ion fraction table: maps sequential index to ionFracs[] index in PhotoIonizationSolver::CellResult.
    // ionFracsIndex = -1 means store n_e (electron density in cm^-3) instead.
    struct IonFracDiagEntry
    {
        int ionFracsIndex;
        const char* name;
    };
    static constexpr IonFracDiagEntry ionFracDiagTable[] = {
        {1, "x_HII"},         // H+ fraction
        {13, "x_NII"},        // N+ fraction
        {20, "x_OI"},         // O0 fraction
        {21, "x_OII"},        // O+ fraction
        {22, "x_OIII"},       // O2+ fraction
        {63, "x_SII"},        // S+ fraction
        {-1, "n_e (cm^-3)"},  // electron density
    };
    static constexpr int numIonFracDiags = sizeof(ionFracDiagTable) / sizeof(ionFracDiagTable[0]);

    // Aggregate ion fraction table: ions tracked as n_HII-weighted global sums.
    // Stored as numbervolumedensity so the framework aggregates sum(x_ion * n_HII * V) across cells.
    // The global mean is then sum(x_ion * n_HII * V) / sum(n_HII * V).
    struct IonFracAggEntry
    {
        int ionFracsIndex;
        const char* name;
    };
    static constexpr IonFracAggEntry ionFracAggTable[] = {
        {13, "nHII_weighted_x_NII"},
        {22, "nHII_weighted_x_OIII"},
        {63, "nHII_weighted_x_SII"},
    };
    static constexpr int numIonFracAggs = sizeof(ionFracAggTable) / sizeof(ionFracAggTable[0]);

    // Density ceiling diagnostic (incremented in initializeSpecificState, logged in isSpecificStateConverged)
    std::atomic<int> densityCeilingCount{0};
    bool densityCeilingLogged{false};

    // PerCell out-of-cross diagnostic (incremented in initializeSpecificState PerCell branch,
    // logged once in isSpecificStateConverged). Cells with imported (n_N, n_C / n_H) implying
    // dN, dC outside the stab cross get endpoint-clamped T/kappa from the stab interpolation.
    std::atomic<int> dNoutOfCrossCount{0};
    std::atomic<int> dCoutOfCrossCount{0};
    std::atomic<int> perCellTotalCount{0};
    bool dCrossLogged{false};

    // ----- Abundance helpers (mass-preserving Gutkin compensation) -----
    // Element-mass-per-H at solar pattern: m_X * n_X/n_H using GASS10 _solarAbundances values.
    // Atomic masses (u) for {C, N, O, Ne, Mg, Si, S, Fe}: {12.011, 14.007, 15.999, 20.180, 24.305, 28.085, 32.065, 55.845}.
    // C=12.011*2.69e-4, N=14.007*6.76e-5, ..., Fe=55.845*3.16e-5. Computed once in C++17 constexpr.
    static constexpr double mCSol = 12.011 * 2.69e-4;
    static constexpr double mNSol = 14.007 * 6.76e-5;
    static constexpr double mOSol = 15.999 * 4.90e-4;
    static constexpr double mNeSol = 20.180 * 8.51e-5;
    static constexpr double mMgSol = 24.305 * 3.98e-5;
    static constexpr double mSiSol = 28.085 * 3.24e-5;
    static constexpr double mSSol = 32.065 * 1.32e-5;
    static constexpr double mFeSol = 55.845 * 3.16e-5;
    static constexpr double mTotSol = mCSol + mNSol + mOSol + mNeSol + mMgSol + mSiSol + mSSol + mFeSol;
    static constexpr double mRestSol = mTotSol - mNSol - mCSol;

    // Gutkin+2016 mass-preserving compensation factor: scales all metals by g so that
    // total Z = solar*Zfrac*mTotSol is preserved when N and C are perturbed by 10^dN, 10^dC.
    // At dN=dC=0, g=1 exactly.
    inline double gutkinCompensation(double dN, double dC)
    {
        return mTotSol / (mNSol * std::pow(10.0, dN) + mCSol * std::pow(10.0, dC) + mRestSol);
    }

    // Gutkin+2016 He enrichment: Y(Z) = 0.2485 + 1.7756 Z; returns n_He/n_H number ratio.
    // Matches the He prescription baked into the multiZ stab generation (see
    // ~/dig_multiZ_stab_pipeline/src/stab_generation/create_stab_compact.py:_he_abundance).
    inline double gutkinHeliumAbundance(double Z)
    {
        const double Y = 0.2485 + 1.7756 * Z;
        const double X = 1.0 - Y - Z;
        return Y / (4.0 * X);
    }

    // Reemission probability indices (following CMACIONIZE nomenclature)
    enum ReemissionChannel {
        Hydrogen = 0,    // Hydrogen Lyman continuum
        HeliumLyC = 1,   // Helium Lyman continuum
        HeliumNpEv = 2,  // He 19.8 eV line
        HeliumTPC = 3,   // Helium two-photon continuum
        HeliumLyA = 4,   // Helium Lyman alpha
    };
    static constexpr int numReemissionChannels = 5;
}

////////////////////////////////////////////////////////////////////

// out-of-line definition of the ReemissionData type forward-declared in the header (see there for
// why): none of its details are needed outside this translation unit
struct DiffuseIonizedGasMix::ReemissionData
{
    bool valid = false;
    double probabilities[numReemissionChannels];
    double cumulativeProbabilities[numReemissionChannels];
    double pHabs;  // Probability of H absorption vs He
};

////////////////////////////////////////////////////////////////////

void DiffuseIonizedGasMix::setupSelfBefore()
{
    EmittingGasMix::setupSelfBefore();

    auto log = find<Log>();

    // Load temperature STAB tables if enabled (multiZ: Standard gains a delta_id axis,
    // Transition stays on the centre slice).
    if (useCloudyTemperature())
    {
        _standardTemperatureTable.open(this, "DiffuseIonizedGas5Bin_Standard_multiZ_Temperature",
                                       "Z(1),n_H(1/m3),logU(1),logR2(1),logR3(1),logR4(1),logR5(1),delta_id(1)",
                                       "logT(K)", true);
        _transitionTemperatureTable.open(this, "DiffuseIonizedGas5Bin_Transition_multiZ_Temperature",
                                         "Z(1),n_H(1/m3),logU(1),logR2(1),logR3(1),logR4(1),logR5(1)", "logT(K)", true);

        // Cache the transition table's upper logU edge, where it meets the standard table.
        Array logUAxis;
        _transitionTemperatureTable.axisArray<2>(logUAxis);
        _transitionLogUMax = logUAxis[logUAxis.size() - 1];
    }

    // Load reemission spectra STAB table if diffuse reemission is enabled
    if (reemissionFraction() > 0.0)
    {
        _reemissionSpectraTable.open(this, "DiffuseIonizedGasReemissionSpectra", "type(1),T(K),lambda(m)", "P(1)",
                                     true);
    }

    // Load opacity STAB tables if enabled (tag-subsampled, K tags in lambda axis;
    // multiZ: Standard carries a delta_id axis, Transition uses only the centre slice).
    if (useCloudyOpacity())
    {
        _standardOpacityTable.open(this, "DiffuseIonizedGas5Bin_Standard_multiZ_Opacity",
                                   "lambda(m),logU(1),logR2(1),logR3(1),logR4(1),logR5(1),Z(1),n_H(1/m3),delta_id(1)",
                                   "kappa(1/m)", true);
        _transitionOpacityTable.open(this, "DiffuseIonizedGas5Bin_Transition_multiZ_Opacity",
                                     "lambda(m),logU(1),logR2(1),logR3(1),logR4(1),logR5(1),Z(1),n_H(1/m3)",
                                     "kappa(1/m)", true);
    }

    // Load DeltaMap table whenever any multiZ stab is in use, and build the dN/dC cross
    // axes for separable bracket+linterp reconstruction. _deltaIdCentre is the delta_id of the
    // (dN=0, dC=0) slice; _dNAxisValues holds the sorted dN grid (with corresponding
    // delta_id indices in _dNAxisDeltaIds), and similarly for dC.
    if (useCloudyTemperature() || useCloudyOpacity())
    {
        _deltaNMapTable.open(this, "DiffuseIonizedGas5Bin_Standard_multiZ_DeltaMap", "delta_id(1)", "delta_N(1)", true);
        _deltaCMapTable.open(this, "DiffuseIonizedGas5Bin_Standard_multiZ_DeltaMap", "delta_id(1)", "delta_C(1)", true);

        const size_t nDelta = _deltaNMapTable.axisSize<0>();
        if (nDelta == 0 || _deltaCMapTable.axisSize<0>() != nDelta)
            throw FATALERROR("DeltaMap table 'DiffuseIonizedGas5Bin_Standard_multiZ_DeltaMap' is empty or has "
                             "mismatched delta_N/delta_C lengths");

        // Treat anything within eps of zero as "exactly on the centre axis". DELTA_SAMPLES
        // values are written exactly, but floating-point round-trips through PTS storage
        // can introduce sub-ULP noise.
        constexpr double eps = 1e-9;
        vector<std::pair<double, int>> dNcandidates;  // dN value, delta_id (entries with dC=0)
        vector<std::pair<double, int>> dCcandidates;  // dC value, delta_id (entries with dN=0)
        for (int deltaId = 0; deltaId < static_cast<int>(nDelta); ++deltaId)
        {
            const double dN = _deltaNMapTable.valueAtIndices(static_cast<size_t>(deltaId));
            const double dC = _deltaCMapTable.valueAtIndices(static_cast<size_t>(deltaId));
            if (std::abs(dN) < eps && std::abs(dC) < eps) _deltaIdCentre = deltaId;
            if (std::abs(dC) < eps) dNcandidates.emplace_back(dN, deltaId);
            if (std::abs(dN) < eps) dCcandidates.emplace_back(dC, deltaId);
        }

        if (_deltaIdCentre < 0)
            throw FATALERROR(
                "DeltaMap table 'DiffuseIonizedGas5Bin_Standard_multiZ_DeltaMap': no centre entry (dN=0, dC=0) found");
        if (dNcandidates.size() < 2)
            throw FATALERROR("DeltaMap table 'DiffuseIonizedGas5Bin_Standard_multiZ_DeltaMap': dN axis has fewer than "
                             "2 grid points");
        if (dCcandidates.size() < 2)
            throw FATALERROR("DeltaMap table 'DiffuseIonizedGas5Bin_Standard_multiZ_DeltaMap': dC axis has fewer than "
                             "2 grid points");

        std::sort(dNcandidates.begin(), dNcandidates.end());
        std::sort(dCcandidates.begin(), dCcandidates.end());
        _dNAxisValues.clear();
        _dNAxisDeltaIds.clear();
        _dCAxisValues.clear();
        _dCAxisDeltaIds.clear();
        for (const auto& p : dNcandidates)
        {
            _dNAxisValues.push_back(p.first);
            _dNAxisDeltaIds.push_back(p.second);
        }
        for (const auto& p : dCcandidates)
        {
            _dCAxisValues.push_back(p.first);
            _dCAxisDeltaIds.push_back(p.second);
        }
        _stabDNmin = _dNAxisValues.front();
        _stabDNmax = _dNAxisValues.back();
        _stabDCmin = _dCAxisValues.front();
        _stabDCmax = _dCAxisValues.back();

        log->info("Stab dN cross: [" + StringUtils::toString(_stabDNmin, 'f', 2) + ", "
                  + StringUtils::toString(_stabDNmax, 'f', 2) + "] dex (" + std::to_string(_dNAxisValues.size())
                  + " points), centre delta_id = " + std::to_string(_deltaIdCentre));
        log->info("Stab dC cross: [" + StringUtils::toString(_stabDCmin, 'f', 2) + ", "
                  + StringUtils::toString(_stabDCmax, 'f', 2) + "] dex (" + std::to_string(_dCAxisValues.size())
                  + " points)");
    }

    _lambdaH = rydbergToWavelength(1.0);   // 1 Ryd = H ionization threshold
    _lambdaHe = rydbergToWavelength(1.8);  // ~1.8 Ryd = He ionization threshold

    // 5-bin wavelength boundaries for ratio calculation
    _lambdaBin1 = rydbergToWavelength(energyBin1Ryd);
    _lambdaBin2 = rydbergToWavelength(energyBin2Ryd);
    _lambdaBin3 = rydbergToWavelength(energyBin3Ryd);
    _lambdaBin4 = rydbergToWavelength(energyBin4Ryd);
    _lambdaBin5 = rydbergToWavelength(energyBin5Ryd);
    _lambdaBin6 = rydbergToWavelength(energyBin6Ryd);

    // _lambdaLow gates opacity (opacityAbs/Sca/Ext + precompute) and selects
    // ionizing photons in calculateIonizationParameter.
    // _lambdaHigh bounds the fallback opacity grid when no RF grid is available.
    _lambdaLow = rydbergToWavelength(1.0) + 1e-10;   // just above 1 Ryd (912 A)
    _lambdaHigh = rydbergToWavelength(6.0) - 2e-10;  // just below 6 Ryd (152 A)

    // Build the gas emission wavelength grid (log-spaced from 900 A): continuum sampling
    // range and the window for extended-inventory lines
    {
        class CustomWavelengthGrid : public DisjointWavelengthGrid
        {
        public:
            void setupSelfBefore() override
            {
                DisjointWavelengthGrid::setupSelfBefore();
                setWavelengthRange(_wl, true);
            }
            void setWavelengths(const Array& wavelengths) { _wl = wavelengths; }

        private:
            Array _wl;
        };

        const int numEmBins = numEmissionWavelengths();
        constexpr double lamMin = 9.0e-8;  // 900 A [m]
        const double lamMax = maxEmissionWavelength();
        Array emLambdav(numEmBins);
        double logMin = std::log10(lamMin);
        double logMax = std::log10(lamMax);
        for (int i = 0; i < numEmBins; i++)
            emLambdav[i] = std::pow(10.0, logMin + (logMax - logMin) * i / (numEmBins - 1));

        auto* grid = new CustomWavelengthGrid();
        grid->setWavelengths(emLambdav);
        const_cast<DiffuseIonizedGasMix*>(this)->addChild(grid);
        grid->setup();
        _emissionWavelengthGrid = grid;
    }

    // Initialize opacity wavelength grid
    initializeOpacityWavelengthGrid();

    if (useCloudyOpacity())
    {
        Array nHRange_opacity;
        _standardOpacityTable.axisArray<7>(nHRange_opacity);
        _nHMinOpacity = nHRange_opacity[0];
    }

    // Initialize emission solver for computing ion fractions from Jv + T.
    // Uses the radiation field wavelength grid for cross-section integration.
    {
        auto config = find<Configuration>();
        auto* rfWavelengthGrid = config->radiationFieldWLG();
        Array rfLambdav, rfDlambdav;
        rfLambdav = rfWavelengthGrid->lambdav();
        rfDlambdav = rfWavelengthGrid->dlambdav();
        _emissionSolver.initialize(rfLambdav, rfDlambdav);

        // Assign each radiation field bin to the energy bin that holds its characteristic wavelength. The stored
        // field per bin is the mean intensity over the bin, so the width-weighted mean over the assigned bins is the
        // exact energy-bin average when the bin borders include the energy-bin boundaries. A radiation field bin
        // that straddles a boundary counts whole in one energy bin; warn so that the grid can be corrected. A border
        // within a relative 1e-4 of a boundary or edge counts as on it (an error below 0.01 % in energy).
        auto straddles = [rfWavelengthGrid](int ell, double lambda) {
            return rfWavelengthGrid->leftBorder(ell) < lambda * (1. - 1e-4)
                   && rfWavelengthGrid->rightBorder(ell) > lambda * (1. + 1e-4);
        };
        auto binText = [rfWavelengthGrid](int ell) {
            return "Radiation field wavelength bin "
                   + StringUtils::toString(rfWavelengthGrid->leftBorder(ell) * 1e6, 'g', 7) + "-"
                   + StringUtils::toString(rfWavelengthGrid->rightBorder(ell) * 1e6, 'g', 7) + " micron";
        };
        const double binLowWl[5] = {_lambdaBin2, _lambdaBin3, _lambdaBin4, _lambdaBin5, _lambdaBin6};
        const double binHighWl[5] = {_lambdaBin1, _lambdaBin2, _lambdaBin3, _lambdaBin4, _lambdaBin5};
        const double boundaries[6] = {_lambdaBin1, _lambdaBin2, _lambdaBin3, _lambdaBin4, _lambdaBin5, _lambdaBin6};
        const double boundariesRyd[6] = {energyBin1Ryd, energyBin2Ryd, energyBin3Ryd,
                                         energyBin4Ryd, energyBin5Ryd, energyBin6Ryd};
        for (int ell = 0; ell < rfWavelengthGrid->numBins(); ell++)
        {
            double lambda = rfWavelengthGrid->wavelength(ell);
            for (int bin = 0; bin < 5; bin++)
            {
                if (lambda >= binLowWl[bin] && lambda <= binHighWl[bin])
                {
                    _energyBinRFBins[bin].push_back(ell);
                    _energyBinRFWidth[bin] += rfWavelengthGrid->effectiveWidth(ell);
                    break;
                }
            }
            for (int k = 0; k < 6; k++)
            {
                if (straddles(ell, boundaries[k]))
                    log->warning(binText(ell) + " straddles the " + StringUtils::toString(boundariesRyd[k], 'g', 3)
                                 + " Ryd boundary of the radiation field ratio bins; place a bin border at "
                                 + StringUtils::toString(boundaries[k] * 1e6, 'g', 7) + " micron");
            }
        }

        // An energy bin without radiation field bins leaves its shape ratio undefined (set to the floor)
        for (int bin = 0; bin < 5; bin++)
        {
            if (_energyBinRFBins[bin].empty())
                log->warning("No radiation field bin has its characteristic wavelength between "
                             + StringUtils::toString(boundariesRyd[bin], 'g', 3) + " and "
                             + StringUtils::toString(boundariesRyd[bin + 1], 'g', 3)
                             + " Ryd; the shape ratio for that range is undefined. The radiation field wavelength "
                               "grid should cover 1 to 6 Ryd ("
                             + StringUtils::toString(_lambdaBin6 * 1e6, 'g', 7) + " to "
                             + StringUtils::toString(_lambdaBin1 * 1e6, 'g', 7) + " micron)");
        }

        // The ion solver evaluates each bin at its characteristic wavelength, so a bin that straddles an ionization
        // edge counts as wholly above or wholly below it. The hydrogen edge matters most: it also decides which bins
        // enter the ionization parameter, and non-ionizing source light in a straddling bin is counted as ionizing.
        const double lambdaHI = Constants::h() * Constants::c() / (VernerCrossSections::HI_eV * Constants::Qelectron());
        const bool sourcesBeyondHI = config->sourceWavelengthRange().max() > lambdaHI;
        string straddledEdges;
        for (int ion = 0; ion < VernerCrossSections::numIons; ion++)
        {
            double E = VernerCrossSections::ionizationPotential(ion);
            double lambdaEdge = Constants::h() * Constants::c() / (E * Constants::Qelectron());
            for (int ell = 0; ell < rfWavelengthGrid->numBins(); ell++)
            {
                if (!straddles(ell, lambdaEdge)) continue;
                if (ion == 0)
                    log->warning(
                        binText(ell) + " straddles the ion solver's hydrogen edge ("
                        + StringUtils::toString(VernerCrossSections::HI_eV, 'g', 4) + " eV); place a bin border at "
                        + StringUtils::toString(lambdaHI * 1e6, 'g', 7) + " micron"
                        + (sourcesBeyondHI ? ". The sources emit beyond this edge, so non-ionizing light in this "
                                             "bin enters the ionization parameter"
                                           : ""));
                else
                    straddledEdges += (straddledEdges.empty() ? "" : ", ") + StringUtils::toString(E, 'g', 4);
            }
        }
        if (!straddledEdges.empty())
            log->warning("Ion solver ionization edges inside a radiation field bin (eV): " + straddledEdges
                         + "; the solver treats each such bin as wholly above or wholly below the edge");
        log->info("DiffuseIonizedGasMix reads the medium number density as total hydrogen plus helium nuclei: "
                  "n_H = n / (1 + y_He)");

        // includeExtendedLines=true requires the full resource set (statistical-equilibrium solver,
        // Case B emissivity tables, extended line registry): a missing resource is a fatal
        // configuration error rather than a silent fallback. includeExtendedLines=false skips all
        // three unconditionally, regardless of whether the resources happen to be installed, so
        // that the fallback lines (legacy P_B recombination, precomputed q_col collisional rates,
        // built-in registry only) stay reachable and testable on their own.
        _gasLineEmission.initialize(this);
        if (includeExtendedLines())
        {
            _gasLineEmission.initializeAtomicModels();
            find<Log>()->info("Collisional line emission uses the statistical-equilibrium solver");

            _gasLineEmission.initializeRecombinationTables();
            find<Log>()->info("Recombination line emission uses the Case B emissivity tables");

            static const char* romans[] = {"I", "II", "III", "IV",   "V",   "VI", "VII", "VIII", "IX",
                                           "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII"};
            static const char* metals[] = {"C", "N", "O", "Ne", "Mg", "Si", "S", "Fe"};
            vector<GasLineEmission::SpeciesSpec> species;
            for (int e = 0; e != 8; ++e)
            {
                int elem = e + 2;
                for (int s = 1; s <= PhotoIonizationSolver::numStages[elem]; ++s)
                    species.push_back(
                        {string(metals[e]) + "_" + romans[s - 1], PhotoIonizationSolver::stageOffset[elem] + s - 1, e});
            }
            _gasLineEmission.initializeExtendedLineRegistry(species);
        }
        else
        {
            find<Log>()->info(
                "Using the built-in fallback lines (legacy P_B recombination, precomputed q_col collisional rates)");
        }

        // select the active lines from the registry (see _activeLines)
        const auto& registry = _gasLineEmission.lineRegistry();
        const Array& emLambdav = _emissionWavelengthGrid->lambdav();
        double emLambdaMin = emLambdav[0];
        double emLambdaMax = emLambdav[emLambdav.size() - 1];
        _activeLines.clear();
        for (int k = 0; k != static_cast<int>(registry.size()); ++k)
        {
            if (k < GasLineEmission::numLines
                || (registry[k].wavelength >= emLambdaMin && registry[k].wavelength <= emLambdaMax))
                _activeLines.push_back(k);
        }
        _numLines = static_cast<int>(_activeLines.size());
        _lineCenters.resize(_numLines);
        _lineMasses.resize(_numLines);
        for (int i = 0; i < _numLines; i++)
        {
            _lineCenters[i] = registry[_activeLines[i]].wavelength;
            _lineMasses[i] = registry[_activeLines[i]].mass;
        }

        // group the collisional lines by atomic model so each species is solved once per cell
        _lineTransitions.assign(_numLines, -1);
        _lineGroup.assign(_numLines, -1);
        for (int i = 0; i < _numLines; ++i)
        {
            int slot = _gasLineEmission.lineModelSlot(_activeLines[i]);
            if (slot < 0) continue;
            _lineTransitions[i] = _gasLineEmission.lineTransition(_activeLines[i]);
            auto it = std::find(_groupSlots.begin(), _groupSlots.end(), slot);
            if (it == _groupSlots.end())
            {
                _groupSlots.push_back(slot);
                _groupLines.emplace_back();
                it = _groupSlots.end() - 1;
            }
            int g = it - _groupSlots.begin();
            _groupLines[g].push_back(i);
            _lineGroup[i] = g;
        }
    }

    // Convert user-facing density ceiling from cm^-3 to m^-3 (0 means disabled)
    if (maxHydrogenDensity() > 0.)
    {
        _maxNH_m3 = maxHydrogenDensity() * 1e6;  // cm^-3 -> m^-3
        log->info("Hydrogen density ceiling: " + StringUtils::toString(maxHydrogenDensity(), 'e', 3) + " cm^-3");
    }
}

////////////////////////////////////////////////////////////////////

bool DiffuseIonizedGasMix::hasExtraSpecificState() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

MaterialMix::DynamicStateType DiffuseIonizedGasMix::hasDynamicMediumState() const
{
    return DynamicStateType::Primary;
}

////////////////////////////////////////////////////////////////////

bool DiffuseIonizedGasMix::hasContinuumEmission() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool DiffuseIonizedGasMix::hasLineEmission() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

Array DiffuseIonizedGasMix::lineEmissionCenters() const
{
    return _lineCenters;
}

////////////////////////////////////////////////////////////////////

Array DiffuseIonizedGasMix::lineEmissionMasses() const
{
    return _lineMasses;
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> DiffuseIonizedGasMix::parameterInfo() const
{
    vector<SnapshotParameter> result;
    if (abundanceMode() == AbundanceMode::PerCell)
    {
        // 8 metal number ratios n_X/n_H, dimensionless. Order MUST match _solarAbundances:
        // C, N, O, Ne, Mg, Si, S, Fe.
        result.push_back(SnapshotParameter::custom("n_C/n_H number ratio"));
        result.push_back(SnapshotParameter::custom("n_N/n_H number ratio"));
        result.push_back(SnapshotParameter::custom("n_O/n_H number ratio"));
        result.push_back(SnapshotParameter::custom("n_Ne/n_H number ratio"));
        result.push_back(SnapshotParameter::custom("n_Mg/n_H number ratio"));
        result.push_back(SnapshotParameter::custom("n_Si/n_H number ratio"));
        result.push_back(SnapshotParameter::custom("n_S/n_H number ratio"));
        result.push_back(SnapshotParameter::custom("n_Fe/n_H number ratio"));
        // 9th column: per-cell helium n_He/n_H.
        result.push_back(SnapshotParameter::custom("n_He/n_H number ratio"));
    }
    return result;
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> DiffuseIonizedGasMix::specificStateVariableInfo() const
{
    // Standard variables first
    vector<StateVariable> result{StateVariable::numberDensity(), StateVariable::metallicity(),
                                 StateVariable::temperature()};

    // Custom variable index counter
    int index = 0;

    // Add custom variable for helium abundance
    const_cast<DiffuseIonizedGasMix*>(this)->_indexHeliumAbundance = index;
    result.push_back(StateVariable::custom(index++, "helium abundance", ""));

    // Add custom variables for neutral fractions
    const_cast<DiffuseIonizedGasMix*>(this)->_indexHNeutralFraction = index;
    result.push_back(StateVariable::custom(index++, "hydrogen neutral fraction", ""));

    const_cast<DiffuseIonizedGasMix*>(this)->_indexHeNeutralFraction = index;
    result.push_back(StateVariable::custom(index++, "helium neutral fraction", ""));

    // Add ionization parameter
    const_cast<DiffuseIonizedGasMix*>(this)->_indexIonizationParameter = index;
    result.push_back(StateVariable::custom(index++, "ionization parameter (logU)", ""));

    // Add 5-bin ratio parameters
    const_cast<DiffuseIonizedGasMix*>(this)->_indexLogR2 = index;
    result.push_back(StateVariable::custom(index++, "log ratio parameter R2", ""));

    const_cast<DiffuseIonizedGasMix*>(this)->_indexLogR3 = index;
    result.push_back(StateVariable::custom(index++, "log ratio parameter R3", ""));

    const_cast<DiffuseIonizedGasMix*>(this)->_indexLogR4 = index;
    result.push_back(StateVariable::custom(index++, "log ratio parameter R4", ""));

    const_cast<DiffuseIonizedGasMix*>(this)->_indexLogR5 = index;
    result.push_back(StateVariable::custom(index++, "log ratio parameter R5", ""));

    // Declared as "numbervolumedensity" so calculateAggregate() automatically computes
    // sum(n_H_ionized * V_cell) across all cells. Neutral cells contribute zero, making
    // the aggregate sensitive only to ionized/front cells (see MediumState.hpp aggregation docs).
    const_cast<DiffuseIonizedGasMix*>(this)->_indexNHIonized = index;
    result.push_back(StateVariable::custom(index++, "ionized hydrogen number density", "numbervolumedensity"));

    // Diagnostic ion fractions from the full PIO solver (always stored)
    const_cast<DiffuseIonizedGasMix*>(this)->_indexFirstIonFractionDiag = index;
    for (int i = 0; i < numIonFracDiags; ++i)
    {
        result.push_back(StateVariable::custom(index++, ionFracDiagTable[i].name, ""));
    }

    // n_HII-weighted ion fraction aggregates for convergence logging
    const_cast<DiffuseIonizedGasMix*>(this)->_indexFirstIonFractionAgg = index;
    for (int i = 0; i < numIonFracAggs; ++i)
    {
        result.push_back(StateVariable::custom(index++, ionFracAggTable[i].name, "numbervolumedensity"));
    }

    // Per-cell metal abundance slots (only in PerCell mode). Order matches _solarAbundances
    // and parameterInfo() snapshot columns: C, N, O, Ne, Mg, Si, S, Fe (number ratios n_X/n_H).
    if (abundanceMode() == AbundanceMode::PerCell)
    {
        const_cast<DiffuseIonizedGasMix*>(this)->_indexFirstMetalAbund = index;
        static const char* names[8] = {"n_C/n_H",  "n_N/n_H",  "n_O/n_H", "n_Ne/n_H",
                                       "n_Mg/n_H", "n_Si/n_H", "n_S/n_H", "n_Fe/n_H"};
        for (int i = 0; i < 8; ++i) result.push_back(StateVariable::custom(index++, names[i], ""));
    }
    else
    {
        const_cast<DiffuseIonizedGasMix*>(this)->_indexFirstMetalAbund = -1;
    }

    const_cast<DiffuseIonizedGasMix*>(this)->_indexFirstOpacityExt = index;
    for (int i = 0; i < _numIonizingOpacityPoints; i++)
    {
        result.push_back(StateVariable::custom(index++, "extinction opacity at wavelength " + std::to_string(i), ""));
    }

    const_cast<DiffuseIonizedGasMix*>(this)->_indexFirstReemissionFraction = index;
    for (int i = 0; i < _numReemissionOpacityPoints; i++)
    {
        result.push_back(StateVariable::custom(
            index++, "effective reemission fraction at wavelength " + std::to_string(_firstReemissionOpacityPoint + i),
            ""));
    }

    return result;
}

////////////////////////////////////////////////////////////////////

// Macro helpers for accessing custom variables in the material state
#define heliumAbundance() custom(_indexHeliumAbundance)
#define hNeutralFraction() custom(_indexHNeutralFraction)
#define heNeutralFraction() custom(_indexHeNeutralFraction)
#define ionizationParameter() custom(_indexIonizationParameter)
#define logR2() custom(_indexLogR2)
#define logR3() custom(_indexLogR3)
#define logR4() custom(_indexLogR4)
#define logR5() custom(_indexLogR5)
#define opacityExtAtIndex(i) custom(_indexFirstOpacityExt + (i))
#define reemissionFractionAtIndex(i) custom(_indexFirstReemissionFraction + (i) - _firstReemissionOpacityPoint)

#define setHeliumAbundance(value) setCustom(_indexHeliumAbundance, (value))
#define setHNeutralFraction(value) setCustom(_indexHNeutralFraction, (value))
#define setHeNeutralFraction(value) setCustom(_indexHeNeutralFraction, (value))
#define setIonizationParameter(value) setCustom(_indexIonizationParameter, (value))
#define setLogR2(value) setCustom(_indexLogR2, (value))
#define setLogR3(value) setCustom(_indexLogR3, (value))
#define setLogR4(value) setCustom(_indexLogR4, (value))
#define setLogR5(value) setCustom(_indexLogR5, (value))
#define nHIonized() custom(_indexNHIonized)
#define setNHIonized(value) setCustom(_indexNHIonized, (value))
#define ionFracDiag(i) custom(_indexFirstIonFractionDiag + (i))
#define setIonFracDiag(i, value) setCustom(_indexFirstIonFractionDiag + (i), (value))
#define ionFracAgg(i) custom(_indexFirstIonFractionAgg + (i))
#define setIonFracAgg(i, value) setCustom(_indexFirstIonFractionAgg + (i), (value))
#define setOpacityExtAtIndex(i, value) setCustom(_indexFirstOpacityExt + (i), (value))
#define setReemissionFractionAtIndex(i, value) \
    setCustom(_indexFirstReemissionFraction + (i) - _firstReemissionOpacityPoint, (value))

////////////////////////////////////////////////////////////////////

void DiffuseIonizedGasMix::initializeSpecificState(MaterialState* state, double metallicity, double temperature,
                                                   const Array& params) const
{
    // If the cell doesn't contain any material, leave all properties at zero
    if (state->numberDensity() > 0.)
    {
        // Hard gate: PerCell metal mode requires a per-cell metallicity column. Falling
        // back to defaultMetallicity() while feeding non-solar N/C/etc. would silently
        // mismatch the stab Z axis vs the analytical abundance vector, so refuse to run.
        if (abundanceMode() == AbundanceMode::PerCell && metallicity < 0.)
            throw FATALERROR("DiffuseIonizedGasMix: PerCell abundance mode requires a per-cell "
                             "metallicity column in the snapshot (map your simulation's gas "
                             "metallicity to the 'Z' parameter); a cell has no imported metallicity.");

        // Resolve cell metallicity (used by all abundance/He modes below)
        const double Z_cell = (metallicity >= 0.) ? metallicity : defaultMetallicity();

        // Default helium: Gutkin Y(Z). PerCell mode overwrites this with the snapshot
        // column below.
        double yHe = gutkinHeliumAbundance(Z_cell);

        // Apply hydrogen density ceiling if configured (clamp total density in-place,
        // so all downstream code automatically uses the capped value)
        if (_maxNH_m3 > 0.)
        {
            double nHm3 = state->numberDensity() / (1.0 + yHe);
            if (nHm3 > _maxNH_m3)
            {
                state->setNumberDensity(_maxNH_m3 * (1.0 + yHe));
                ++densityCeilingCount;
            }
        }

        // Set metallicity and temperature
        state->setMetallicity(Z_cell);
        state->setTemperature(temperature >= 0. ? temperature : defaultTemperature());

        // Per-cell metal abundances (PerCell mode only). Snapshot params come in the order
        // declared by parameterInfo(): C, N, O, Ne, Mg, Si, S, Fe, [He].
        int paramIdx = 0;
        if (abundanceMode() == AbundanceMode::PerCell)
        {
            for (int i = 0; i < 8; ++i) state->setCustom(_indexFirstMetalAbund + i, params[paramIdx++]);

            // Diagnostic: track cells whose imported (n_N, n_C / n_H) imply dN, dC outside
            // the stab cross. Outside-cross cells get endpoint-clamped T/kappa silently.
            // Counts are logged once at the first convergence check (see isSpecificStateConverged).
            const double eps = 1e-30;
            const double Zfrac = Z_cell / _solarZ;
            const double n_C_imported = state->custom(_indexFirstMetalAbund + 0);
            const double n_N_imported = state->custom(_indexFirstMetalAbund + 1);
            const double dN_cell = std::log10(std::max(n_N_imported, eps) / std::max(_solarAbundances[1] * Zfrac, eps));
            const double dC_cell = std::log10(std::max(n_C_imported, eps) / std::max(_solarAbundances[0] * Zfrac, eps));
            if (dN_cell < _stabDNmin || dN_cell > _stabDNmax) ++dNoutOfCrossCount;
            if (dC_cell < _stabDCmin || dC_cell > _stabDCmax) ++dCoutOfCrossCount;
            ++perCellTotalCount;
        }

        // PerCell mode reads the 9th snapshot column (n_He/n_H) here, overwriting the
        // Gutkin Y(Z) baseline resolved above. SolarScaled keeps the Gutkin value.
        if (abundanceMode() == AbundanceMode::PerCell) yHe = params[paramIdx++];
        state->setHeliumAbundance(yHe);

        // Initial neutral fractions - start with quasi-ionized gas
        state->setHNeutralFraction(minNeutralFraction);   // Essentially zero
        state->setHeNeutralFraction(minNeutralFraction);  // Essentially zero

        // Initialize ionized H density reflecting the fully-ionized initial condition,
        // so the first previousAggregate correctly represents the starting state.
        double n_H_init = convertTotalDensityToHydrogenDensity(state->numberDensity(), yHe);
        state->setNHIonized(n_H_init * (1.0 - minNeutralFraction));

        // Initialize radiation field parameters to defaults
        state->setIonizationParameter(_logFloor);
        state->setLogR2(_logFloor);
        state->setLogR3(_logFloor);
        state->setLogR4(_logFloor);
        state->setLogR5(_logFloor);

        // Initialize diagnostic ion fractions and aggregates to zero
        for (int i = 0; i < numIonFracDiags; ++i) state->setIonFracDiag(i, 0.0);
        for (int i = 0; i < numIonFracAggs; ++i) state->setIonFracAgg(i, 0.0);

        // Initialize opacity arrays with zeros
        for (int i = 0; i < _numIonizingOpacityPoints; i++) state->setOpacityExtAtIndex(i, 0.0);
        for (int i = _firstReemissionOpacityPoint; i < _firstReemissionOpacityPoint + _numReemissionOpacityPoints; i++)
            state->setReemissionFractionAtIndex(i, 0.0);
    }
}

////////////////////////////////////////////////////////////////////

UpdateStatus DiffuseIonizedGasMix::updateSpecificState(MaterialState* state, const Array& Jv) const
{
    UpdateStatus status;

    // Only update if the cell contains material
    if (state->numberDensity() > 0.)
    {
        // For convergence check: save previous state
        double oldT = state->temperature();
        double oldLogU = state->ionizationParameter();

        // Calculate and store ionization parameter
        double n_total = state->numberDensity();
        double He_abundance = state->heliumAbundance();
        double n_H = convertTotalDensityToHydrogenDensity(n_total, He_abundance);
        double newLogU = calculateIonizationParameter(Jv, n_H);
        state->setIonizationParameter(newLogU);

        // Calculate and store 5-bin ratio parameters
        double logR2_val, logR3_val, logR4_val, logR5_val;
        calculate5BinRatioParameters(Jv, logR2_val, logR3_val, logR4_val, logR5_val);
        state->setLogR2(logR2_val);
        state->setLogR3(logR3_val);
        state->setLogR4(logR4_val);
        state->setLogR5(logR5_val);

        // Update temperature from STAB before ionization balance (PIO solver needs T)
        if (useCloudyTemperature())
        {
            updateTemperatureFromStab(state);
        }

        // Full ionization balance via PIO solver (replaces standalone H/He solver).
        // Uses converged Jv and current T to compute all ion fractions self-consistently.
        constexpr double cm3PerM3 = 1e6;
        double nH_cgs = n_H / cm3PerM3;  // m^-3 -> cm^-3
        double T = state->temperature();
        double solverIonFracs[PhotoIonizationSolver::totalStages] = {};
        bool hasSolverIonFracs = false;

        if (T > 0.)
        {
            double metalAbundances[8];
            buildPerCellAbundances(state, metalAbundances);

            auto result =
                _emissionSolver.solveIonizationAtFixedT(Jv, nH_cgs, He_abundance, metalAbundances, T, nullptr);
            for (int s = 0; s < PhotoIonizationSolver::totalStages; ++s) solverIonFracs[s] = result.ionFracs[s];
            hasSolverIonFracs = true;

            // H/He neutral fractions from the full solver
            double h0 = std::max(minNeutralFraction, std::min(result.ionFracs[0], maxNeutralFraction));
            double he0 = std::max(minNeutralFraction, std::min(result.ionFracs[2], maxNeutralFraction));
            state->setHNeutralFraction(h0);
            state->setHeNeutralFraction(he0);

            // Store diagnostic ion fractions
            for (int i = 0; i < numIonFracDiags; ++i)
            {
                if (ionFracDiagTable[i].ionFracsIndex >= 0)
                    state->setIonFracDiag(i, result.ionFracs[ionFracDiagTable[i].ionFracsIndex]);
                else
                    state->setIonFracDiag(i, result.ne);
            }

            // Store n_HII-weighted ion fraction aggregates: x_ion * n_HII [m^-3]
            double nHII = n_H * (1.0 - h0);
            for (int i = 0; i < numIonFracAggs; ++i)
                state->setIonFracAgg(i, result.ionFracs[ionFracAggTable[i].ionFracsIndex] * nHII);
        }
        else
        {
            state->setHNeutralFraction(1.0);
            state->setHeNeutralFraction(1.0);
            for (int i = 0; i < numIonFracDiags; ++i) state->setIonFracDiag(i, 0.0);
            for (int i = 0; i < numIonFracAggs; ++i) state->setIonFracAgg(i, 0.0);
        }

        // Pre-compute opacity arrays
        precomputeOpacityArrays(state, Jv, hasSolverIonFracs ? solverIonFracs : nullptr);

        // Update ionized H density for global convergence criterion
        state->setNHIonized(n_H * (1.0 - state->hNeutralFraction()));

        // Check per-cell convergence: cell fails if either T or logU changed significantly.
        // Only check logU for cells above the emission floor; cells below have noisy logU
        // from MC noise but contribute zero emission, so their convergence is irrelevant.
        double newT = state->temperature();
        double changeT = (oldT > 0.) ? std::abs(newT - oldT) / oldT : 0.;
        bool logURelevant = (newLogU > _minLogUEmit);
        double changeLogU =
            logURelevant ? std::abs(newLogU - oldLogU) / std::max(std::abs(oldLogU), minDivisionSafe) : 0.;

        if (changeT > maxChangeInIonization() || changeLogU > maxChangeInIonization())
            status.updateNotConverged();
        else
            status.updateConverged();
    }

    return status;
}

////////////////////////////////////////////////////////////////////

bool DiffuseIonizedGasMix::isSpecificStateConverged(int /*numCells*/, int numUpdated, int numNotConverged,
                                                    MaterialState* currentAggregate,
                                                    MaterialState* previousAggregate) const
{
    // Calculate fraction of converged cells (0-1)
    // Use numUpdated instead of numCells to only count cells that contain material; when this
    // component currently has no active cells at all, there is nothing that could have failed to
    // converge, so treat it as vacuously converged instead of dividing by zero
    double fractionNotConverged =
        numUpdated > 0 ? static_cast<double>(numNotConverged) / static_cast<double>(numUpdated) : 0.;
    double convergedFraction = 1.0 - fractionNotConverged;

    // Check standard convergence criterion
    bool standardConverged = fractionNotConverged <= maxFractionNotConvergedCells();

    // Check stability-based convergence: converged fraction unchanged for 3 consecutive iterations
    bool stabilityConverged = false;
    _convergedFractionHistory.push_back(convergedFraction);

    // Keep only the last N iterations
    if (_convergedFractionHistory.size() > _convergenceHistorySize)
    {
        _convergedFractionHistory.erase(_convergedFractionHistory.begin());
    }

    // Check that there is enough history and that the fraction is stable.
    if (_convergedFractionHistory.size() >= _convergenceHistorySize)
    {
        bool isStable = true;
        for (size_t i = 1; i < _convergedFractionHistory.size(); ++i)
        {
            double change = std::abs(_convergedFractionHistory[i] - _convergedFractionHistory[i - 1]);
            if (change > stabilityConvergenceThreshold())
            {
                isStable = false;
                break;
            }
        }
        stabilityConverged = isStable;
    }

    // Global criterion: relative change in total ionized H between consecutive iterations.
    // Neutral cells contribute zero to the aggregate, so this is sensitive only to ionized/front
    // cells. Prevents fake-convergence when trivially-neutral bulk cells dominate the per-cell count.
    double currentTotalIonizedH = currentAggregate->nHIonized();
    double previousTotalIonizedH = previousAggregate->nHIonized();
    double globalChange =
        std::abs(currentTotalIonizedH - previousTotalIonizedH) / std::max(currentTotalIonizedH, minDivisionSafe);
    bool globalConverged = globalChange <= maxChangeInGlobalIonizedH();

    // Overall convergence: (per-cell OR plateau) AND global ionization structure stable.
    // The global gate prevents the OR from being exploited by trivially-neutral cells.
    bool converged = (standardConverged || stabilityConverged) && globalConverged;

    auto log = find<Log>();

    // Log density ceiling count once (cells were capped in initializeSpecificState)
    if (!densityCeilingLogged && _maxNH_m3 > 0.)
    {
        int nCapped = densityCeilingCount.load();
        if (nCapped > 0)
            log->info("Density ceiling: " + std::to_string(nCapped)
                      + " cells capped to nH = " + StringUtils::toString(maxHydrogenDensity(), 'e', 3) + " cm^-3");
        densityCeilingLogged = true;
    }

    // Log PerCell out-of-cross summary once (cells flagged in initializeSpecificState).
    if (!dCrossLogged && abundanceMode() == AbundanceMode::PerCell)
    {
        int total = perCellTotalCount.load();
        if (total > 0)
        {
            int nOutN = dNoutOfCrossCount.load();
            int nOutC = dCoutOfCrossCount.load();
            log->info("PerCell out-of-cross: dN " + std::to_string(nOutN) + "/" + std::to_string(total) + " ("
                      + StringUtils::toString(100.0 * nOutN / total, 'f', 2) + "%);" + " dC " + std::to_string(nOutC)
                      + "/" + std::to_string(total) + " (" + StringUtils::toString(100.0 * nOutC / total, 'f', 2) + "%)"
                      + "  [outside-cross cells get endpoint-clamped T/kappa]");
        }
        dCrossLogged = true;
    }

    // Compact single-line summary: header + 3 criteria flags.
    // Plateau history is appended when at least 2 samples are available.
    string plateauHistory;
    if (_convergedFractionHistory.size() >= 2)
    {
        plateauHistory = " ";
        for (size_t i = 0; i < _convergedFractionHistory.size(); ++i)
        {
            plateauHistory += StringUtils::toString(_convergedFractionHistory[i] * 100., 'f', 1) + "%";
            if (i < _convergedFractionHistory.size() - 1) plateauHistory += "->";
        }
    }
    log->info("DiffuseIonizedGasMix convergence: " + string(converged ? "CONVERGED" : "NOT CONVERGED") + " | per-cell "
              + StringUtils::toString(convergedFraction * 100., 'f', 1) + "%/"
              + StringUtils::toString((1.0 - maxFractionNotConvergedCells()) * 100., 'f', 1) + "% "
              + string(standardConverged ? "PASS" : "FAIL") + " | plateau" + plateauHistory + " "
              + string(stabilityConverged ? "PASS" : "FAIL") + " | global dnHII "
              + StringUtils::toString(globalChange * 100., 'f', 2) + "%/"
              + StringUtils::toString(maxChangeInGlobalIonizedH() * 100., 'f', 2) + "% "
              + string(globalConverged ? "PASS" : "FAIL"));

    // Log n_HII-weighted mean ion fractions: <x_ion> = sum(x_ion * n_HII * V) / sum(n_HII * V)
    if (currentTotalIonizedH > 0.)
    {
        string ionStr = "  Ion fractions (n_HII-weighted):";
        for (int i = 0; i < numIonFracAggs; ++i)
        {
            double currentVal = currentAggregate->ionFracAgg(i);
            double meanIonFrac = currentVal / currentTotalIonizedH;
            ionStr += " " + string(ionFracAggTable[i].name).substr(16)  // strip "nHII_weighted_" prefix
                      + "=" + StringUtils::toString(meanIonFrac, 'e', 3);
            // Also log iteration-to-iteration change
            double previousVal = previousAggregate->ionFracAgg(i);
            if (previousVal > 0.)
            {
                double prevMean = previousVal / std::max(previousTotalIonizedH, minDivisionSafe);
                double relChange = std::abs(meanIonFrac - prevMean) / std::max(prevMean, minDivisionSafe);
                ionStr += " (d=" + StringUtils::toString(relChange * 100., 'f', 1) + "%)";
            }
        }
        log->info(ionStr);
    }

    // Reset history if converged (for next round of iterations)
    if (converged)
    {
        _convergedFractionHistory.clear();
    }

    return converged;
}

////////////////////////////////////////////////////////////////////

double DiffuseIonizedGasMix::mass() const
{
    // Indicative mean mass per particle at the default metallicity; uses Gutkin Y(Z)
    // to stay consistent with the stab and analytical-solver baseline.
    const double yHe = gutkinHeliumAbundance(defaultMetallicity());
    return Constants::Mproton() * (1.0 + 4.0 * yHe) / (1.0 + yHe);
}

////////////////////////////////////////////////////////////////////

double DiffuseIonizedGasMix::sectionAbs(double /*lambda*/) const
{
    return 0.0;
}

////////////////////////////////////////////////////////////////////

bool DiffuseIonizedGasMix::hasScattering() const
{
    return reemissionFraction() > 0.0;
}

////////////////////////////////////////////////////////////////////

bool DiffuseIonizedGasMix::hasScatteringDispersion() const
{
    return reemissionFraction() > 0.0;
}

////////////////////////////////////////////////////////////////////

double DiffuseIonizedGasMix::sectionSca(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double DiffuseIonizedGasMix::sectionExt(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double DiffuseIonizedGasMix::opacityAbs(double lambda, const MaterialState* state, const PhotonPacket*) const
{
    // Return zero opacity if no material or outside ionizing wavelength range
    if (state->numberDensity() <= 0. || lambda > _lambdaLow) return 0.;

    // Use pre-computed opacity array with interpolation
    MaterialState* nonConstState = const_cast<MaterialState*>(state);

    // Use interpolation with absorption opacity type (0)
    return interpolateOpacityFromState(lambda, nonConstState, 0);
}

////////////////////////////////////////////////////////////////////

double DiffuseIonizedGasMix::opacitySca(double lambda, const MaterialState* state, const PhotonPacket*) const
{
    // Scattering opacity equals total opacity minus absorption opacity when reemission is enabled
    if (reemissionFraction() <= 0.0 || state->numberDensity() <= 0.) return 0.;

    // Reemission only in 1-6 Ryd; precomputed scattering is zero beyond
    if (lambda > _lambdaLow) return 0.;

    // Use pre-computed opacity array with interpolation
    MaterialState* nonConstState = const_cast<MaterialState*>(state);

    // Use interpolation with scattering opacity type (1)
    return interpolateOpacityFromState(lambda, nonConstState, 1);
}

////////////////////////////////////////////////////////////////////

double DiffuseIonizedGasMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket*) const
{
    // Return zero opacity if no material or outside ionizing wavelength range
    if (state->numberDensity() <= 0. || lambda > _lambdaLow) return 0.;

    // Use pre-computed opacity array with interpolation
    MaterialState* nonConstState = const_cast<MaterialState*>(state);

    // Use interpolation with extinction opacity type (2)
    return interpolateOpacityFromState(lambda, nonConstState, 2);
}

////////////////////////////////////////////////////////////////////

bool DiffuseIonizedGasMix::peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda,
                                             Direction bfkobs, Direction /*bfky*/, const MaterialState* state,
                                             const PhotonPacket* pp) const
{
    // Initialize to no scattering
    I = Q = U = V = 0.;

    if (reemissionFraction() <= 0.0 || state->numberDensity() <= 0.) return false;

    // Draw random reemission channel unless already stored
    auto scatinfo = const_cast<PhotonPacket*>(pp)->getScatteringInfo();
    if (!scatinfo->valid)
    {
        int channel = selectReemissionChannel(state, lambda, pp);
        scatinfo->valid = true;
        scatinfo->species = channel;

        // Determine new wavelength based on channel
        double newLambda = lambda;

        switch (channel)
        {
            case ReemissionChannel::Hydrogen: newLambda = sampleHydrogenLymanContinuum(state->temperature()); break;
            case ReemissionChannel::HeliumLyC: newLambda = sampleHeliumLymanContinuum(state->temperature()); break;
            case ReemissionChannel::HeliumNpEv: newLambda = Constants::c() / helium19p8EvFreq; break;
            case ReemissionChannel::HeliumTPC: newLambda = sampleHeliumTwoPhotonContinuum(state->temperature()); break;
            case ReemissionChannel::HeliumLyA:
                // This should not happen as selectReemissionChannel turns this into absorbed or TPC
                scatinfo->lambda = 0.;
                return false;
            default:
                // No reemission, photon absorbed
                scatinfo->lambda = 0.;
                return false;
        }

        scatinfo->lambda = newLambda;
    }

    // For reemission, emission is isotropic and unpolarized
    if (scatinfo->lambda > 0.)
    {
        lambda = scatinfo->lambda;
        I = 1.0;         // Bias weight is 1 for isotropic emission
        Q = U = V = 0.;  // Unpolarized

        // Apply Doppler shift due to bulk velocity if present
        if (state->bulkVelocity().norm() > 0.)
        {
            lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfkobs, state->bulkVelocity());
        }
    }
    return false;
}

////////////////////////////////////////////////////////////////////

void DiffuseIonizedGasMix::performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const
{
    if (reemissionFraction() <= 0.0 || state->numberDensity() <= 0.) return;

    // Get or calculate scattering info
    auto scatinfo = pp->getScatteringInfo();
    if (!scatinfo->valid)
    {
        int channel = selectReemissionChannel(state, lambda, pp);
        scatinfo->valid = true;
        scatinfo->species = channel;

        // Determine new wavelength based on channel
        double newLambda = lambda;

        switch (channel)
        {
            case ReemissionChannel::Hydrogen: newLambda = sampleHydrogenLymanContinuum(state->temperature()); break;
            case ReemissionChannel::HeliumLyC: newLambda = sampleHeliumLymanContinuum(state->temperature()); break;
            case ReemissionChannel::HeliumNpEv: newLambda = Constants::c() / helium19p8EvFreq; break;
            case ReemissionChannel::HeliumTPC: newLambda = sampleHeliumTwoPhotonContinuum(state->temperature()); break;
            case ReemissionChannel::HeliumLyA:
                // This should not happen as selectReemissionChannel turns this into absorbed or TPC
                newLambda = 0.;
                break;
            default:
                // No reemission, photon absorbed
                newLambda = 0.;
                break;
        }

        scatinfo->lambda = newLambda;
    }

    // Execute the scattering/reemission event
    if (scatinfo->lambda > 0.)
    {
        // Reemission: isotropic direction, new wavelength
        Direction newDirection = random()->direction();
        pp->setUnpolarized();  // Reemission is unpolarized
        pp->scatter(newDirection, state->bulkVelocity(), scatinfo->lambda);
    }
    else
    {
        // Absorption: terminate photon
    }
}

////////////////////////////////////////////////////////////////////

DisjointWavelengthGrid* DiffuseIonizedGasMix::emissionWavelengthGrid() const
{
    return _emissionWavelengthGrid;
}

////////////////////////////////////////////////////////////////////

DiffuseIonizedGasMix::TableSelection DiffuseIonizedGasMix::selectTable(double logU) const
{
    // The transition table's upper logU edge was cached in setupSelfBefore.
    const double blendWidth = transitionBlendWidth();

    // The tables share data only at their common edge, _transitionLogUMax, so blend
    // there and nowhere else. Below it Transition applies and clamps at its own
    // floor; Standard would clamp to the shared edge, handing near-neutral gas the
    // opacity of gas at the ionised-neutral boundary. The opacity path is not gated
    // by _minLogUEmit.
    if (logU <= _transitionLogUMax - blendWidth) return TableSelection::Transition;
    if (logU > _transitionLogUMax + blendWidth) return TableSelection::Standard;
    return TableSelection::Blend;
}

////////////////////////////////////////////////////////////////////

double DiffuseIonizedGasMix::getBlendingWeight(double logU) const
{
    // Returns weight for transition table (0 = standard only, 1 = transition only).
    // The transition table's upper logU edge was cached in setupSelfBefore.
    const double blendWidth = transitionBlendWidth();

    // Blend only at the tables' shared edge; see selectTable().
    if (logU > _transitionLogUMax - blendWidth && logU <= _transitionLogUMax + blendWidth)
    {
        double weight = ((_transitionLogUMax + blendWidth) - logU) / (2.0 * blendWidth);
        return std::max(0.0, std::min(1.0, weight));
    }
    if (logU <= _transitionLogUMax - blendWidth) return 1.0;
    return 0.0;
}

////////////////////////////////////////////////////////////////////

Array DiffuseIonizedGasMix::emissionSpectrum(const MaterialState* state, const Array& Jv) const
{
    constexpr double cm3PerM3 = 1e6;

    int numWavelengths = _emissionWavelengthGrid->numBins();
    Array emission(numWavelengths);

    if (state->numberDensity() <= 0.) return emission;
    if (state->ionizationParameter() <= _minLogUEmit) return emission;

    // Cell properties in CGS
    double T = state->temperature();
    if (T <= 0.) return emission;
    double yHe = state->heliumAbundance();
    double nHm3 = state->numberDensity() / (1.0 + yHe);
    double nH = nHm3 / cm3PerM3;  // [cm^-3]
    double V_cm3 = state->volume() * 1e6;

    // Solve ion fractions from converged Jv + T
    double metalAbundances[8];
    buildPerCellAbundances(state, metalAbundances);

    auto result = _emissionSolver.solveIonizationAtFixedT(Jv, nH, yHe, metalAbundances, T, nullptr);

    double ne = result.ne;
    if (ne <= 0.) return emission;

    // Ion number densities for continuum [cm^-3]
    double nHII = nH * result.ionFracs[1];                                                  // HII
    double nHeII = nH * yHe * result.ionFracs[PhotoIonizationSolver::stageOffset[1] + 1];   // HeII
    double nHeIII = nH * yHe * result.ionFracs[PhotoIonizationSolver::stageOffset[1] + 2];  // HeIII

    // Continuum emission (free-bound + free-free + two-photon)
    const Array& lambdav = _emissionWavelengthGrid->lambdav();
    for (int i = 0; i < numWavelengths; i++)
    {
        emission[i] = GasContinuumEmission::continuumLuminosity(lambdav[i], T, ne, nHII, nHeII, nHeIII, V_cm3);
    }

    return emission;
}

////////////////////////////////////////////////////////////////////

Array DiffuseIonizedGasMix::lineEmissionSpectrum(const MaterialState* state, const Array& Jv) const
{
    constexpr double cm3PerM3 = 1e6;

    Array luminosities(_numLines);

    if (state->numberDensity() <= 0.) return luminosities;
    if (state->ionizationParameter() <= _minLogUEmit) return luminosities;

    // Cell properties in CGS
    double T = state->temperature();
    if (T <= 0.) return luminosities;
    double yHe = state->heliumAbundance();
    double nHm3 = state->numberDensity() / (1.0 + yHe);
    double nH = nHm3 / cm3PerM3;  // [cm^-3]
    double V_cm3 = state->volume() * 1e6;

    // Metal abundances per abundanceMode (SolarScaled: solar pattern * Zfrac * g; PerCell: snapshot)
    double metalAbundances[8];
    buildPerCellAbundances(state, metalAbundances);

    // Solve ion fractions from converged Jv + T
    auto result = _emissionSolver.solveIonizationAtFixedT(Jv, nH, yHe, metalAbundances, T, nullptr);

    double ne = result.ne;
    if (ne <= 0.) return luminosities;

    double nHI = nH * result.ionFracs[0];

    // recombining ion densities [cm^-3] for the recombination lines
    double nHII = nH * result.ionFracs[1];
    double nHeII = nH * yHe * result.ionFracs[PhotoIonizationSolver::stageOffset[1] + 1];
    double nHeIII = nH * yHe * result.ionFracs[PhotoIonizationSolver::stageOffset[1] + 2];

    // Gamma_HI for the legacy P_B path, from the same Jv the solver was called with so that it
    // is consistent with the ion fractions just computed
    double gamma[PhotoIonizationSolver::totalStages];
    _emissionSolver.computePhotoionizationRates(Jv, gamma);
    double gammaHI = gamma[0];
    const auto& registry = _gasLineEmission.lineRegistry();

    // collisional lines served by an atomic model: one level-population solve per species,
    // in the nebular limit (electron collisions at T and ne, no radiative pumping)
    for (size_t g = 0; g != _groupSlots.size(); ++g)
    {
        const auto& line = registry[_activeLines[_groupLines[g].front()]];
        double nIon = nH * metalAbundances[line.elementIndex] * result.ionFracs[line.carrierIonIndex];
        GasLineEmission::Environment env;
        env.Tkin = T;
        env.nTotal = nIon * 1e6;    // cm^-3 -> m^-3
        env.nPartner = {ne * 1e6};  // cm^-3 -> m^-3
        auto pops = _gasLineEmission.solveLevelPopulations(_gasLineEmission.atomicModel(_groupSlots[g]), env);
        auto eps = _gasLineEmission.lineEmissivities(_gasLineEmission.atomicModel(_groupSlots[g]), pops);
        for (int i : _groupLines[g]) luminosities[i] = eps[_lineTransitions[i]] * V_cm3 * 1e-6;  // W m^-3 x m^3
    }

    // remaining lines: recombination lines (Case B tables, legacy P_B fallback for H) and
    // collisional lines on the legacy table path
    for (int i = 0; i < _numLines; ++i)
    {
        if (_lineGroup[i] >= 0) continue;
        int k = _activeLines[i];
        const auto& line = registry[k];
        if (line.carrierIonIndex >= 0)
        {
            double nIon = nH * metalAbundances[line.elementIndex] * result.ionFracs[line.carrierIonIndex];
            luminosities[i] = _gasLineEmission.collisionalLineLuminosity(k, T, ne, nIon, V_cm3);
        }
        else
        {
            double nIon = line.isHeIRecomb ? nHeII : (line.isHeIIRecomb ? nHeIII : nHII);
            luminosities[i] = _gasLineEmission.recombinationLineLuminosity(k, T, ne, nIon, gammaHI, nHI, V_cm3);
        }
    }

    // optional per-cell relative floor on the extended-inventory lines
    if (lineLuminosityFloor() > 0.)
    {
        double Lmax = 0.;
        for (int i = 0; i < _numLines; ++i) Lmax = max(Lmax, luminosities[i]);
        double floor = lineLuminosityFloor() * Lmax;
        for (int i = 0; i < _numLines; ++i)
            if (_activeLines[i] >= GasLineEmission::numLines && luminosities[i] < floor) luminosities[i] = 0.;
    }

    return luminosities;
}

////////////////////////////////////////////////////////////////////

double DiffuseIonizedGasMix::indicativeTemperature(const MaterialState* state, const Array&) const
{
    return state->temperature();
}

////////////////////////////////////////////////////////////////////

void DiffuseIonizedGasMix::updateTemperatureFromStab(MaterialState* state) const
{
    // Get current parameters
    double Z = state->metallicity();
    double n_total = state->numberDensity();
    double He_abundance = state->heliumAbundance();
    double n_H = convertTotalDensityToHydrogenDensity(n_total, He_abundance);
    double logU = state->ionizationParameter();
    double logR2 = state->logR2();
    double logR3 = state->logR3();
    double logR4 = state->logR4();
    double logR5 = state->logR5();

    // If all R values are invalid (at initial values) or ionization parameter below emission floor,
    // set temperature to floor. emissionSpectrum() also gates on _minLogUEmit, keeping state consistent.
    if ((logR2 <= _logFloor && logR3 <= _logFloor && logR4 <= _logFloor && logR5 <= _logFloor) || logU <= _minLogUEmit)
    {
        state->setTemperature(_minTEmit);
        return;
    }

    // Determine which table(s) to use
    TableSelection choice = selectTable(logU);
    double logT = 0.0;

    // Standard temperature is reconstructed via bracket+linterp+separable on the dN/dC cross:
    //   logT(dN, dC) = logT_dN(dN) + logT_dC(dC) - logT_centre.
    // The Transition table has no delta_id axis (uses centre slice by construction).
    const double centreDeltaId = static_cast<double>(_deltaIdCentre);
    double dN, dC;
    computeCellDeltas(state, dN, dC);
    int dNlo, dNhi, dClo, dChi;
    double wN, wC;
    bracketDeltaAxis(dN, _dNAxisValues, _dNAxisDeltaIds, dNlo, dNhi, wN);
    bracketDeltaAxis(dC, _dCAxisValues, _dCAxisDeltaIds, dClo, dChi, wC);

    auto stdTempSeparable = [&]() -> double {
        const double logT_dN_lo =
            _standardTemperatureTable(Z, n_H, logU, logR2, logR3, logR4, logR5, static_cast<double>(dNlo));
        const double logT_dN_hi =
            _standardTemperatureTable(Z, n_H, logU, logR2, logR3, logR4, logR5, static_cast<double>(dNhi));
        const double logT_dC_lo =
            _standardTemperatureTable(Z, n_H, logU, logR2, logR3, logR4, logR5, static_cast<double>(dClo));
        const double logT_dC_hi =
            _standardTemperatureTable(Z, n_H, logU, logR2, logR3, logR4, logR5, static_cast<double>(dChi));
        const double logT_ctr = _standardTemperatureTable(Z, n_H, logU, logR2, logR3, logR4, logR5, centreDeltaId);
        const double logT_dN = (1.0 - wN) * logT_dN_lo + wN * logT_dN_hi;
        const double logT_dC = (1.0 - wC) * logT_dC_lo + wC * logT_dC_hi;
        return logT_dN + logT_dC - logT_ctr;
    };

    switch (choice)
    {
        case TableSelection::Standard: logT = stdTempSeparable(); break;

        case TableSelection::Transition:
            logT = _transitionTemperatureTable(Z, n_H, logU, logR2, logR3, logR4, logR5);
            break;

        case TableSelection::Blend:
        {
            double w = getBlendingWeight(logU);
            double logT_std = stdTempSeparable();
            double logT_trans = _transitionTemperatureTable(Z, n_H, logU, logR2, logR3, logR4, logR5);

            // Linear blend in log-space (temperature varies smoothly in log)
            logT = (1.0 - w) * logT_std + w * logT_trans;
            break;
        }
    }

    double newTemperature = std::pow(10.0, logT);
    state->setTemperature(newTemperature);
}

////////////////////////////////////////////////////////////////////

double DiffuseIonizedGasMix::calculateIonizationParameter(const Array& Jv, double nH) const
{
    /**
     * Calculate ionization parameter directly from radiation field:
     * U = (ionizing photon flux) / (n_H * c)
     *
     * Since we use step functions with energy conservation by construction,
     * we can calculate U directly without Cloudy field reconstruction.
     */

    constexpr double h = Constants::h();
    constexpr double c = Constants::c();

    // Get radiation field wavelength grid
    auto config = find<Configuration>();
    auto rfwlg = config->radiationFieldWLG();

    // Ionizing photon flux: phi = sum over bins with characteristic wavelength above 1 Ryd of
    // 4pi * J_lambda * lambda * dlambda / (h*c). Jv holds the mean intensity over each bin, so the sum is exact
    // when a bin border sits at the ionization edge.
    double ionizing_photon_flux = 0.;
    int numIonizingBins = 0;
    for (int i = 0; i < rfwlg->numBins(); i++)
    {
        double lambda = rfwlg->wavelength(i);
        if (lambda <= _lambdaLow)
        {
            ionizing_photon_flux += 4. * M_PI * Jv[i] * lambda * rfwlg->effectiveWidth(i) / (h * c);
            numIonizingBins++;
        }
    }

    // If no ionizing radiation, return the sentinel value
    if (numIonizingBins == 0) return _logFloor;

    // Calculate ionization parameter: U = phi / (n_H * c)
    if (nH < 1e-20) nH = 1e-20;  // Prevent division by zero

    double U = ionizing_photon_flux / (nH * c);
    double logU = log10(std::max(U, 1e-99));

    // Apply reasonable bounds
    const double maxLogU = 5.0;  // Very high ionization
    logU = std::max(_logFloor, std::min(logU, maxLogU));

    return logU;
}

////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
// Helper functions for the diffuse reemision

void DiffuseIonizedGasMix::interpolateReemissionSpectrum(int spectrumType, double temperature,
                                                         vector<double>& wavelengths,
                                                         vector<double>& cumulativeDist) const
{
    // Get wavelength grid from STAB table
    Array lambdaArray;
    _reemissionSpectraTable.axisArray<2>(lambdaArray);

    wavelengths.resize(lambdaArray.size());
    cumulativeDist.resize(lambdaArray.size());

    for (size_t i = 0; i < lambdaArray.size(); ++i)
    {
        wavelengths[i] = lambdaArray[i];
    }

    // Interpolate cumulative distribution for each wavelength bin
    for (size_t i = 0; i < wavelengths.size(); ++i)
    {
        // Interpolate: table(spectrum_type, temperature, wavelength) -> cumulative_probability
        cumulativeDist[i] = _reemissionSpectraTable(static_cast<double>(spectrumType), temperature, wavelengths[i]);
    }

    // Ensure monotonicity and proper normalization
    for (size_t i = 1; i < cumulativeDist.size(); ++i)
    {
        if (cumulativeDist[i] < cumulativeDist[i - 1])
        {
            cumulativeDist[i] = cumulativeDist[i - 1];
        }
    }

    // Normalize to ensure final value is 1.0
    if (cumulativeDist.back() > 0.0)
    {
        double normalization = 1.0 / cumulativeDist.back();
        for (double& value : cumulativeDist)
        {
            value *= normalization;
        }
    }
}

////////////////////////////////////////////////////////////////////

double DiffuseIonizedGasMix::sampleFromTemperatureDependentSpectrum(int spectrumType, double temperature,
                                                                    Random* random) const
{
    // Cache the last (spectrumType, temperature) CDF so consecutive reemission events at
    // the same temperature reuse it instead of rebuilding it from STAB lookups.
    thread_local int cachedSpectrumType = -1;
    thread_local double cachedTemperature = 0.;
    thread_local vector<double> cachedWavelengths;
    thread_local vector<double> cachedCumulativeDist;

    if (spectrumType != cachedSpectrumType || temperature != cachedTemperature)
    {
        interpolateReemissionSpectrum(spectrumType, temperature, cachedWavelengths, cachedCumulativeDist);
        cachedSpectrumType = spectrumType;
        cachedTemperature = temperature;
    }
    const vector<double>& wavelengths = cachedWavelengths;
    const vector<double>& cumulativeDist = cachedCumulativeDist;

    // Sample using inverse transform method
    double x = random->uniform();

    // Find position in cumulative distribution
    auto it = std::lower_bound(cumulativeDist.begin(), cumulativeDist.end(), x);

    if (it == cumulativeDist.end())
    {
        return wavelengths.back();
    }

    size_t index = std::distance(cumulativeDist.begin(), it);

    if (index == 0)
    {
        return wavelengths[0];
    }

    // Linear interpolation between adjacent bins
    double lambda1 = wavelengths[index - 1];
    double lambda2 = wavelengths[index];
    double prob1 = cumulativeDist[index - 1];
    double prob2 = cumulativeDist[index];

    if (prob2 != prob1)
    {
        double fraction = (x - prob1) / (prob2 - prob1);
        return lambda1 + fraction * (lambda2 - lambda1);
    }
    else
    {
        return lambda1;
    }
}

////////////////////////////////////////////////////////////////////

void DiffuseIonizedGasMix::calculateReemissionProbabilities(const MaterialState* state, double lambda,
                                                            ReemissionData& data) const
{
    data.valid = false;

    if (state->numberDensity() <= 0.)
    {
        return;
    }

    // Following CMACIONIZE PhysicalDiffuseReemissionHandler::set_reemission_probabilities
    // Many of the following are in CGS, but ratios are used eventually

    const double T = state->temperature();
    const double T4 = T * 1e-4;  // Temperature in units of 10^4 K

    // Hydrogen reemission probability
    // Wood, Mathis & Ercolano (2004), equation (24)
    const double alpha_1_H = 1.58e-13 * std::pow(T4, -0.53);
    const double alpha_total = 4.18e-13 * std::pow(T4, -0.7);  // value used in CMACIONIZE
    const double prob_H = alpha_1_H / alpha_total;

    // Helium reemission probabilities
    // Wood, Mathis & Ercolano (2004), equation (25)
    // The published values differ from the values used in CMACIONIZE;
    // the CMACIONIZE values are adopted here.
    const double alpha_1_He = 1.54e-13 * std::pow(T4, -0.486);   // He Lyman continuum
    const double alpha_e_2tS = 2.1e-13 * std::pow(T4, -0.381);   // 23S
    const double alpha_e_2sS = 2.06e-14 * std::pow(T4, -0.451);  // 21S
    const double alpha_e_2sP = 4.17e-14 * std::pow(T4, -0.695);  // 21P

    const double alphaHe_total = alpha_1_He + alpha_e_2tS + alpha_e_2sS + alpha_e_2sP;

    // Store individual probabilities (normalized, assume that there are only these channels)
    data.probabilities[ReemissionChannel::Hydrogen] = prob_H;
    data.probabilities[ReemissionChannel::HeliumLyC] = alpha_1_He / alphaHe_total;
    data.probabilities[ReemissionChannel::HeliumNpEv] = alpha_e_2tS / alphaHe_total;
    data.probabilities[ReemissionChannel::HeliumTPC] = alpha_e_2sS / alphaHe_total;
    data.probabilities[ReemissionChannel::HeliumLyA] = alpha_e_2sP / alphaHe_total;

    // Create cumulative probabilities for helium channels
    data.cumulativeProbabilities[ReemissionChannel::HeliumLyC] = data.probabilities[ReemissionChannel::HeliumLyC];
    data.cumulativeProbabilities[ReemissionChannel::HeliumNpEv] =
        data.cumulativeProbabilities[ReemissionChannel::HeliumLyC] + data.probabilities[ReemissionChannel::HeliumNpEv];
    data.cumulativeProbabilities[ReemissionChannel::HeliumTPC] =
        data.cumulativeProbabilities[ReemissionChannel::HeliumNpEv] + data.probabilities[ReemissionChannel::HeliumTPC];
    data.cumulativeProbabilities[ReemissionChannel::HeliumLyA] = 1.0;  // Ensures total probability = 1

    // Calculate H vs He absorption probability using wavelength-dependent cross-sections
    const double h0 = state->hNeutralFraction();
    const double he0 = state->heNeutralFraction();
    const double AHe = state->heliumAbundance();

    if (h0 > 0. || (he0 > 0. && AHe > 0.))
    {
        // wavelength-dependent cross-sections
        const double frequency = Constants::c() / lambda;
        const double sigmaH = getHydrogenCrossSection(frequency);

        // Helium cross-section only if wavelength can ionize helium
        double sigmaHe = 0.0;
        if (lambda <= _lambdaHe)  // Only if energetic enough to ionize helium
        {
            sigmaHe = getHeliumCrossSection(frequency);
        }

        // Calculate absorption probabilities weighted by neutral fractions and abundances
        const double nH0_sigma = h0 * sigmaH;
        const double nHe0_sigma = AHe * he0 * sigmaHe;

        const double total_sigma = nH0_sigma + nHe0_sigma;
        if (total_sigma > 0.)
        {
            data.pHabs = nH0_sigma / total_sigma;
        }
        else
        {
            // This can happen when wavelength is outside ionization range due to slight mismatches
            // Just mark as invalid and return
            data.valid = false;
            return;
        }
    }
    else
    {
        throw FATALERROR("DIG: Diffuse reemission called for fully ionized gas with no neutral atoms. "
                         "This should not happen as scattering opacity should be zero. ");
    }

    data.valid = true;
}

////////////////////////////////////////////////////////////////////

int DiffuseIonizedGasMix::selectReemissionChannel(const MaterialState* state, double lambda,
                                                  const PhotonPacket* /*pp*/) const
{
    if (state->numberDensity() <= 0.)
    {
        return -1;
    }

    ReemissionData data;
    calculateReemissionProbabilities(state, lambda, data);

    if (!data.valid)
    {
        return -1;  // Absorption (no reemission)
    }

    // First decide if absorbed by H or He
    double x = random()->uniform();

    if (x <= data.pHabs)
    {
        // Absorbed by hydrogen
        x = random()->uniform();
        if (x <= data.probabilities[ReemissionChannel::Hydrogen])
        {
            return ReemissionChannel::Hydrogen;
        }
        else
        {
            return -1;  // Absorbed, no reemission
        }
    }
    else
    {
        // Absorbed by helium - select reemission channel
        x = random()->uniform();

        for (int i = ReemissionChannel::HeliumLyC; i <= ReemissionChannel::HeliumLyA; ++i)
        {
            if (x <= data.cumulativeProbabilities[i])
            {
                // Special handling for two-photon continuum
                if (i == ReemissionChannel::HeliumTPC)
                {
                    // Only 56% of two-photon events produce H-ionizing radiation
                    x = random()->uniform();
                    if (x < heliumTpcHIonizingFraction)
                    {
                        return ReemissionChannel::HeliumTPC;
                    }
                    else
                    {
                        return -1;  // Absorbed
                    }
                }

                // Special handling for Helium Lyman alpha (on-the-spot absorption)
                if (i == ReemissionChannel::HeliumLyA)
                {
                    const double h0 = state->hNeutralFraction();
                    const double he0 = state->heNeutralFraction();
                    const double T = state->temperature();

                    // Calculate on-the-spot absorption probability
                    const double sqrtT_nH0 = std::sqrt(T) * h0;
                    const double pHots = sqrtT_nH0 / (sqrtT_nH0 + heliumLyaOtsCoeff * he0);

                    x = random()->uniform();
                    if (x < pHots)
                    {
                        // Absorbed on-the-spot by hydrogen
                        x = random()->uniform();
                        if (x <= data.probabilities[ReemissionChannel::Hydrogen])
                        {
                            return ReemissionChannel::Hydrogen;
                        }
                        else
                        {
                            return -1;  // Absorbed
                        }
                    }
                    else
                    {
                        // Converted to helium two-photon continuum
                        x = random()->uniform();
                        if (x < heliumTpcHIonizingFraction)
                        {
                            return ReemissionChannel::HeliumTPC;
                        }
                        else
                        {
                            return -1;  // Absorbed
                        }
                    }
                }

                return i;
            }
        }

        return -1;  // Should not reach here, but return absorption as fallback
    }
}

////////////////////////////////////////////////////////////////////

double DiffuseIonizedGasMix::sampleHydrogenLymanContinuum(double temperature) const
{
    // from temperature-dependent STAB table
    const int hLymanSpectrumType = 0;
    return sampleFromTemperatureDependentSpectrum(hLymanSpectrumType, temperature, random());
}

////////////////////////////////////////////////////////////////////

double DiffuseIonizedGasMix::sampleHeliumLymanContinuum(double temperature) const
{
    // from temperature-dependent STAB table
    const int heLymanSpectrumType = 1;
    return sampleFromTemperatureDependentSpectrum(heLymanSpectrumType, temperature, random());
}

////////////////////////////////////////////////////////////////////

double DiffuseIonizedGasMix::sampleHeliumTwoPhotonContinuum(double temperature) const
{
    // from temperature-dependent STAB table
    const int heTwoPhotonSpectrumType = 2;
    return sampleFromTemperatureDependentSpectrum(heTwoPhotonSpectrumType, temperature, random());
}

////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////

void DiffuseIonizedGasMix::calculate5BinRatioParameters(const Array& Jv, double& logR2, double& logR3, double& logR4,
                                                        double& logR5) const
{
    // 5-bin ratio parameters from the radiation field: R_i = J_i / J_1 for i = 2..5,
    // with J_i the average intensity in bin i

    // Initialize ratios to the sentinel value; STAB will clamp valid ratios into bounds.
    logR2 = _logFloor;
    logR3 = _logFloor;
    logR4 = _logFloor;
    logR5 = _logFloor;

    // Calculate bin averages
    double binAverages[5];
    calculateBinAverages(Jv, binAverages);

    // Calculate energy content in each bin for validation
    double binWidths[5] = {
        energyBin2Ryd - energyBin1Ryd,  // Bin 1: 1.00 - 1.80 Ryd
        energyBin3Ryd - energyBin2Ryd,  // Bin 2: 1.80 - 2.58 Ryd
        energyBin4Ryd - energyBin3Ryd,  // Bin 3: 2.58 - 3.52 Ryd
        energyBin5Ryd - energyBin4Ryd,  // Bin 4: 3.52 - 4.00 Ryd
        energyBin6Ryd - energyBin5Ryd   // Bin 5: 4.00 - 6.00 Ryd
    };

    double totalEnergy = 0.0;
    double binEnergies[5];
    for (int i = 0; i < 5; i++)
    {
        binEnergies[i] = binAverages[i] * binWidths[i];
        totalEnergy += binEnergies[i];
    }

    // Define threshold for significant energy fraction
    const double minEnergyFraction = 1e-32;
    const double minSignificantIntensity = 1e-99;

    // Check if bins are significant
    bool isSignificantBin[5];
    for (int i = 0; i < 5; i++)
    {
        double energyFraction = (totalEnergy > 0.0) ? (binEnergies[i] / totalEnergy) : 0.0;
        isSignificantBin[i] = (energyFraction >= minEnergyFraction) && (binAverages[i] > minSignificantIntensity);
    }

    // Calculate ratios only if bins are significant
    if (isSignificantBin[0] && binAverages[0] > minSignificantIntensity)
    {
        // R2 = <J2>/<J1>
        if (isSignificantBin[1])
        {
            double R2 = binAverages[1] / binAverages[0];
            R2 = std::max(1e-99, R2);  // Ensure positive
            logR2 = log10(R2);
        }

        // R3 = <J3>/<J1>
        if (isSignificantBin[2])
        {
            double R3 = binAverages[2] / binAverages[0];
            R3 = std::max(1e-99, R3);  // Ensure positive
            logR3 = log10(R3);
        }

        // R4 = <J4>/<J1>
        if (isSignificantBin[3])
        {
            double R4 = binAverages[3] / binAverages[0];
            R4 = std::max(1e-99, R4);  // Ensure positive
            logR4 = log10(R4);
        }

        // R5 = <J5>/<J1>
        if (isSignificantBin[4])
        {
            double R5 = binAverages[4] / binAverages[0];
            R5 = std::max(1e-99, R5);  // Ensure positive
            logR5 = log10(R5);
        }
    }
}

////////////////////////////////////////////////////////////////////

void DiffuseIonizedGasMix::calculateBinAverages(const Array& Jv, double* binAverages) const
{
    /**
     * Calculate bin averages for the 5-bin method.
     * Bin 1: 1.0 - 1.80 Ryd
     * Bin 2: 1.8 - 2.58 Ryd
     * Bin 3: 2.58 - 3.52 Ryd
     * Bin 4: 3.52 - 4.00 Ryd
     * Bin 5: 4.00 - 6.00 Ryd
     */

    // Jv holds the mean intensity over each radiation field bin, so the average over an energy bin is the
    // width-weighted mean over the radiation field bins assigned to it during setup.
    // An energy bin without radiation field bins gets the minimal value.
    auto rfwlg = find<Configuration>()->radiationFieldWLG();
    for (int bin = 0; bin < 5; bin++)
    {
        double sum = 0.;
        for (int ell : _energyBinRFBins[bin]) sum += Jv[ell] * rfwlg->effectiveWidth(ell);
        double average = _energyBinRFWidth[bin] > 0. ? sum / _energyBinRFWidth[bin] : 0.;
        binAverages[bin] = std::max(average, 1e-99);
    }
}

////////////////////////////////////////////////////////////////////

double DiffuseIonizedGasMix::rydbergToWavelength(double energy_ryd) const
{
    constexpr double h = Constants::h();
    constexpr double c = Constants::c();
    constexpr double rydberg = rydbergEnergyEv * Constants::Qelectron();

    return h * c / (energy_ryd * rydberg);
}

////////////////////////////////////////////////////////////////////

double DiffuseIonizedGasMix::getHydrogenCrossSection(double frequency) const
{
    // convert frequency to photon energy in eV, then use the canonical Verner+96 fit
    const double E_eV = Constants::h() * frequency / Constants::Qelectron();
    return VernerCrossSections::sigmaHI(E_eV) * 1e-4;  // cm^2 -> m^2
}

////////////////////////////////////////////////////////////////////

double DiffuseIonizedGasMix::getHeliumCrossSection(double frequency) const
{
    // convert frequency to photon energy in eV, then use the canonical Verner+96 fit
    const double E_eV = Constants::h() * frequency / Constants::Qelectron();
    return VernerCrossSections::sigmaHeI(E_eV) * 1e-4;  // cm^2 -> m^2
}

////////////////////////////////////////////////////////////////////

double DiffuseIonizedGasMix::interpolateOpacityFromState(double lambda, MaterialState* state, int opacityType) const
{
    const int numWavelengths = _opacityWavelengthGrid.size();
    if (lambda < _opacityWavelengthGrid[0] || lambda > _opacityWavelengthGrid[numWavelengths - 1]) return 0.;
    if (opacityType < 0 || opacityType > 2) return 0.;

    // stored opacity of the requested type at a given grid index; the absorption and scattering parts are rebuilt
    // from the extinction opacity and the effective re-emission fraction with the same arithmetic that split them
    auto stored = [this, state, opacityType](int k) {
        if (k >= _numIonizingOpacityPoints) return 0.;
        double ext = state->opacityExtAtIndex(k);
        if (opacityType == 2) return ext;
        bool inReemission =
            k >= _firstReemissionOpacityPoint && k < _firstReemissionOpacityPoint + _numReemissionOpacityPoints;
        if (opacityType == 0) return inReemission ? ext * (1.0 - state->reemissionFractionAtIndex(k)) : ext;
        return inReemission ? ext * state->reemissionFractionAtIndex(k) : 0.;
    };

    // only the two grid points bracketing lambda are read, with the same bracketing as NR::clampedValue
    int i = NR::locate(_opacityWavelengthGrid, lambda);
    if (i < 0) return stored(0);
    if (i >= numWavelengths - 1) return stored(numWavelengths - 1);
    double x1 = _opacityWavelengthGrid[i];
    double x2 = _opacityWavelengthGrid[i + 1];
    double y1 = stored(i);
    double y2 = stored(i + 1);

    // grid points beyond the ionization edge hold zero opacity by construction; a log-log interpolation towards
    // such a point returns zero, which would let ionizing photons between the last ionizing grid point and the
    // H I edge (13.598434 eV, NIST) cross the gas unabsorbed, so these photons get the opacity of the last
    // ionizing grid point instead
    constexpr double lambdaEdgeHI = Constants::h() * Constants::c() / (13.598434 * Constants::Qelectron());
    if (y2 <= 0. && y1 > 0. && x2 > _lambdaLow && lambda <= lambdaEdgeHI) return y1;

    return NR::interpolateLogLog(lambda, x1, x2, y1, y2);
}

////////////////////////////////////////////////////////////////////

void DiffuseIonizedGasMix::initializeOpacityWavelengthGrid() const
{
    // Only initialize once
    if (_opacityWavelengthGrid.size() > 0) return;

    // Get the radiation field wavelength grid
    auto config = find<Configuration>();
    auto rfwlg = config->radiationFieldWLG();

    if (rfwlg)
    {
        // Use the radiation field wavelength grid
        _opacityWavelengthGrid = rfwlg->lambdav();
    }
    else
    {
        // Fallback: create a simple wavelength grid in the relevant range
        const int numPoints = 100;
        _opacityWavelengthGrid.resize(numPoints);

        // Log-space grid from _lambdaHigh to _lambdaLow
        double logLambdaMin = std::log10(_lambdaHigh);
        double logLambdaMax = std::log10(_lambdaLow);
        double deltaLogLambda = (logLambdaMax - logLambdaMin) / (numPoints - 1);

        for (int i = 0; i < numPoints; i++)
        {
            _opacityWavelengthGrid[i] = std::pow(10.0, logLambdaMin + i * deltaLogLambda);
        }

        auto log = find<Log>();
        log->warning("DiffuseIonizedGasMix: No radiation field grid available, using fallback opacity grid ("
                     + std::to_string(numPoints) + " points)");
    }

    // Grid points that can hold nonzero opacity (ionizing) and those in the re-emission range; the grid is sorted
    // by increasing wavelength, so each set is a contiguous range of indices
    const int numWavelengths = _opacityWavelengthGrid.size();
    _numIonizingOpacityPoints = 0;
    while (_numIonizingOpacityPoints < numWavelengths
           && _opacityWavelengthGrid[_numIonizingOpacityPoints] <= _lambdaLow)
        _numIonizingOpacityPoints++;
    _firstReemissionOpacityPoint = 0;
    _numReemissionOpacityPoints = 0;
    if (reemissionFraction() > 0.0)
    {
        while (_firstReemissionOpacityPoint < numWavelengths
               && _opacityWavelengthGrid[_firstReemissionOpacityPoint] < _lambdaBin5)
            _firstReemissionOpacityPoint++;
        while (_firstReemissionOpacityPoint + _numReemissionOpacityPoints < numWavelengths
               && _opacityWavelengthGrid[_firstReemissionOpacityPoint + _numReemissionOpacityPoints] <= _lambdaH)
            _numReemissionOpacityPoints++;
    }
    find<Log>()->info("DiffuseIonizedGasMix: opacities stored at " + std::to_string(_numIonizingOpacityPoints)
                      + " ionizing grid points, re-emission fractions at " + std::to_string(_numReemissionOpacityPoints)
                      + " points (" + std::to_string(numWavelengths) + " radiation field grid points)");
}

////////////////////////////////////////////////////////////////////

void DiffuseIonizedGasMix::precomputeOpacityArrays(MaterialState* state, const Array& Jv,
                                                   const double* solverIonFracs) const
{
    // grid points beyond the ionizing range hold zero opacity and are not stored
    const int numWavelengths = _numIonizingOpacityPoints;

    // Get actual density
    double n_total = state->numberDensity();
    double He_abundance = state->heliumAbundance();
    double n_H = convertTotalDensityToHydrogenDensity(n_total, He_abundance);
    double logU = state->ionizationParameter();

    // Scale down opacity for cells below the table's minimum density
    // to avoid an over-opaque medium at low densities (StoredTable clamps to boundary)
    // At higher densities, the opacity table is held flat for consistency with the emission tables (see the class docstring).
    double densityScalingFactor = 1.0;
    if (useCloudyOpacity() && n_H < _nHMinOpacity)
    {
        densityScalingFactor = n_H / _nHMinOpacity;
    }

    // Below the logU range of the Cloudy opacity tables (the emission floor, where the temperature is set to the
    // floor value and the lines are switched off), the table lookup would clamp to the edge of the table, which
    // describes partly ionized gas. These cells use the analytical opacity from the ion fractions of the
    // ionization-balance solver instead, as do all cells when the Cloudy opacity is not used.
    const bool useAnalyticalOpacity = !useCloudyOpacity() || logU <= _minLogUEmit;

    // Ion fractions for the analytical opacity (used across all wavelengths): taken from the ionization-balance
    // solution of this state update when available, otherwise solved here
    double analyticalIonFracs[PhotoIonizationSolver::totalStages] = {};
    double analyticalAbundances[8] = {};
    double nH_cgs_opacity = 0.;
    bool analyticalOpacityReady = false;
    if (useAnalyticalOpacity)
    {
        double T = state->temperature();
        double yHe = state->heliumAbundance();
        nH_cgs_opacity = n_H / 1e6;  // m^-3 to cm^-3
        // Per-cell abundances per abundanceMode (SolarScaled: solar * Zfrac * g; PerCell: snapshot).
        buildPerCellAbundances(state, analyticalAbundances);

        if (solverIonFracs)
        {
            for (int s = 0; s < PhotoIonizationSolver::totalStages; ++s) analyticalIonFracs[s] = solverIonFracs[s];
            analyticalOpacityReady = true;
        }
        else if (T > 0. && logU > -98.)
        {
            auto result =
                _emissionSolver.solveIonizationAtFixedT(Jv, nH_cgs_opacity, yHe, analyticalAbundances, T, nullptr);
            for (int s = 0; s < PhotoIonizationSolver::totalStages; ++s) analyticalIonFracs[s] = result.ionFracs[s];
            analyticalOpacityReady = true;
        }
    }

    // dN/dC bracket+weights for the Standard opacity reconstruction. Cell-level (independent
    // of wavelength), so cache once outside the loop. Only used when the Standard table is
    // queried (TableSelection::Standard or Blend); the Transition table has no delta_id axis.
    const double centreDeltaId = static_cast<double>(_deltaIdCentre);
    int dNlo = _deltaIdCentre, dNhi = _deltaIdCentre;
    int dClo = _deltaIdCentre, dChi = _deltaIdCentre;
    double wN = 0., wC = 0.;
    if (!useAnalyticalOpacity)
    {
        double dN, dC;
        computeCellDeltas(state, dN, dC);
        bracketDeltaAxis(dN, _dNAxisValues, _dNAxisDeltaIds, dNlo, dNhi, wN);
        bracketDeltaAxis(dC, _dCAxisValues, _dCAxisDeltaIds, dClo, dChi, wC);
    }
    const double dNloD = static_cast<double>(dNlo), dNhiD = static_cast<double>(dNhi);
    const double dCloD = static_cast<double>(dClo), dChiD = static_cast<double>(dChi);

    for (int i = 0; i < numWavelengths; i++)
    {
        double lambda = _opacityWavelengthGrid[i];
        double totalOpacity = 0.0;

        // Calculate reemission probability once (used by both STAB and non-STAB paths)
        double probReemission = 0.0;
        bool inReemissionRange = (lambda >= _lambdaBin5 && lambda <= _lambdaH);

        if (inReemissionRange)
        {
            // Get reemission data (includes all probabilities)
            ReemissionData data;
            calculateReemissionProbabilities(state, lambda, data);

            // if data is invalid, leave probReemission at its 0 default
            if (data.valid)
            {
                // Calculate hydrogen reemission probability
                double hydrogenScatProb = data.probabilities[ReemissionChannel::Hydrogen];

                // Calculate helium reemission probability (if wavelength can ionize helium)
                double heliumScatProb = 0.0;
                if (lambda <= _lambdaHe)
                {
                    heliumScatProb += data.probabilities[ReemissionChannel::HeliumLyC];
                    heliumScatProb += data.probabilities[ReemissionChannel::HeliumNpEv];
                    heliumScatProb += data.probabilities[ReemissionChannel::HeliumTPC] * heliumTpcHIonizingFraction;

                    // Handle Helium Lyman alpha on-the-spot absorption
                    const double h0 = state->hNeutralFraction();
                    const double he0 = state->heNeutralFraction();
                    const double T = state->temperature();
                    const double sqrtT_nH0 = std::sqrt(T) * h0;
                    const double pHots = sqrtT_nH0 / (sqrtT_nH0 + heliumLyaOtsCoeff * he0);
                    heliumScatProb += data.probabilities[ReemissionChannel::HeliumLyA] * (1.0 - pHots);
                }

                // Weighted average reemission probability based on which species absorbs
                probReemission = data.pHabs * hydrogenScatProb + (1.0 - data.pHabs) * heliumScatProb;
            }
        }

        if (!useAnalyticalOpacity)
        {
            if (lambda <= _lambdaLow)
            {
                double opacityValue = 0.0;
                double logR2 = state->logR2();
                double logR3 = state->logR3();
                double logR4 = state->logR4();
                double logR5 = state->logR5();
                double Z = state->metallicity();

                TableSelection choice = selectTable(logU);
                // Standard opacity reconstructed via bracket+linterp+separable on the dN/dC
                // cross. Linear blend in log-kappa to keep the (smooth in log) opacity behaviour:
                //   log10(kappa)(dN, dC) = log10(kappa_dN(dN)) + log10(kappa_dC(dC)) - log10(kappa_centre).
                auto stdOpacitySeparable = [&]() -> double {
                    const double k_dN_lo =
                        _standardOpacityTable(lambda, logU, logR2, logR3, logR4, logR5, Z, n_H, dNloD);
                    const double k_dN_hi =
                        _standardOpacityTable(lambda, logU, logR2, logR3, logR4, logR5, Z, n_H, dNhiD);
                    const double k_dC_lo =
                        _standardOpacityTable(lambda, logU, logR2, logR3, logR4, logR5, Z, n_H, dCloD);
                    const double k_dC_hi =
                        _standardOpacityTable(lambda, logU, logR2, logR3, logR4, logR5, Z, n_H, dChiD);
                    const double k_ctr =
                        _standardOpacityTable(lambda, logU, logR2, logR3, logR4, logR5, Z, n_H, centreDeltaId);
                    const double eps = 1e-99;
                    const double logK_dN =
                        (1.0 - wN) * std::log10(std::max(k_dN_lo, eps)) + wN * std::log10(std::max(k_dN_hi, eps));
                    const double logK_dC =
                        (1.0 - wC) * std::log10(std::max(k_dC_lo, eps)) + wC * std::log10(std::max(k_dC_hi, eps));
                    const double logK_ctr = std::log10(std::max(k_ctr, eps));
                    return std::pow(10.0, logK_dN + logK_dC - logK_ctr);
                };

                switch (choice)
                {
                    case TableSelection::Standard: opacityValue = stdOpacitySeparable(); break;
                    case TableSelection::Transition:
                        opacityValue = _transitionOpacityTable(lambda, logU, logR2, logR3, logR4, logR5, Z, n_H);
                        break;
                    case TableSelection::Blend:
                    {
                        double w = getBlendingWeight(logU);
                        double stdVal = stdOpacitySeparable();
                        double transVal = _transitionOpacityTable(lambda, logU, logR2, logR3, logR4, logR5, Z, n_H);
                        double logStd = log10(std::max(stdVal, 1e-99));
                        double logTrans = log10(std::max(transVal, 1e-99));
                        opacityValue = pow(10.0, (1.0 - w) * logStd + w * logTrans);
                        break;
                    }
                }

                totalOpacity = opacityValue;

                // Apply density scaling if outside table bounds (linear scaling)
                totalOpacity *= densityScalingFactor;
            }
        }
        else  // analytical opacity from the ion fractions of the ionization-balance solver
        {
            if (lambda <= _lambdaLow && analyticalOpacityReady)
            {
                // opacityAbs returns cm^-1; convert to m^-1
                totalOpacity = _emissionSolver.opacityAbs(lambda, analyticalIonFracs, nH_cgs_opacity, He_abundance,
                                                          analyticalAbundances)
                               * 100.0;
            }
        }

        // the absorption part is totalOpacity * (1 - f) and the scattering part totalOpacity * f, with f the
        // effective re-emission fraction (interpolateOpacityFromState)
        state->setOpacityExtAtIndex(i, totalOpacity);
        if (i >= _firstReemissionOpacityPoint && i < _firstReemissionOpacityPoint + _numReemissionOpacityPoints)
            state->setReemissionFractionAtIndex(i, probReemission * reemissionFraction());
    }
}

////////////////////////////////////////////////////////////////////

double DiffuseIonizedGasMix::convertTotalDensityToHydrogenDensity(double n_total, double He_abundance) const
{
    // For H+He mix: n_total = n_H * (1 + He/H)
    double conversion_factor = 1.0 + He_abundance;
    return n_total / conversion_factor;
}

////////////////////////////////////////////////////////////////////

void DiffuseIonizedGasMix::buildPerCellAbundances(const MaterialState* state, double abund[8]) const
{
    if (abundanceMode() == AbundanceMode::PerCell)
    {
        // Read per-cell metal abundances from snapshot-filled state slots.
        // Order matches _solarAbundances and parameterInfo(): C, N, O, Ne, Mg, Si, S, Fe.
        for (int i = 0; i < 8; ++i) abund[i] = state->custom(_indexFirstMetalAbund + i);
    }
    else
    {
        // SolarScaled mode: mass-preserving Gutkin compensation with N and C scalings.
        // linscale[i] carries the per-element linear scaling on top of solar*Zfrac*g.
        // Order MUST match _solarAbundances: {C, N, O, Ne, Mg, Si, S, Fe}; only C and N
        // are user-tunable, the rest follow the solar pattern (factor 1).
        const double Zfrac = state->metallicity() / _solarZ;
        const double dN = std::log10(solarScalingN());
        const double dC = std::log10(solarScalingC());
        const double scaledZfrac = Zfrac * gutkinCompensation(dN, dC);
        const double linscale[8] = {solarScalingC(), solarScalingN(), 1, 1, 1, 1, 1, 1};
        for (int i = 0; i < 8; ++i) abund[i] = _solarAbundances[i] * scaledZfrac * linscale[i];
    }

    // Apply the Gunasekera (2023) gas-phase depletion factors to bring the
    // analytical ionization and line-emission step into agreement with the
    // composition Cloudy used to generate the temperature and opacity tables.
    // See class docstring (section "Abundances and gas-phase depletion") for
    // the conventions and the source of the per-element factors.
    for (int i = 0; i < 8; ++i) abund[i] *= _gunasekeraDepletion[i];
}

////////////////////////////////////////////////////////////////////

void DiffuseIonizedGasMix::computeCellDeltas(const MaterialState* state, double& dN, double& dC) const
{
    if (abundanceMode() == AbundanceMode::PerCell)
    {
        // Recover the Gutkin-recipe (dN, dC) that, paired with the imported Z and the
        // mass-preserving compensation g(dN, dC), reproduces the imported (n_C, n_N).
        // The stab axis labels recipe inputs, not physical (n_X / solar*Zfrac) ratios:
        // Cloudy at axis point (dN, dC) saw n_N/n_H = solar_N * Zfrac * g * 10^dN.
        // Reading the ratio directly biases the stab lookup by log10(g) at non-zero dN
        // or dC, so invert the recipe instead.
        //
        // Defining a = n_N_imp / (solar_N*Zfrac), b = n_C_imp / (solar_C*Zfrac), the
        // recipe constraints a = g*10^dN, b = g*10^dC together with the definition of g
        // give the closed form g = (M_tot - M_N*a - M_C*b) / M_rest. Exact on recipe
        // cells; for arbitrary cells (MEGATRON), projects onto the (dN, dC) that would
        // reproduce the imported N and C abundances under the recipe family.
        constexpr double eps = 1e-30;
        const double Zfrac = state->metallicity() / _solarZ;
        const double n_C_imported = state->custom(_indexFirstMetalAbund + 0);
        const double n_N_imported = state->custom(_indexFirstMetalAbund + 1);
        const double a = std::max(n_N_imported, eps) / std::max(_solarAbundances[1] * Zfrac, eps);
        const double b = std::max(n_C_imported, eps) / std::max(_solarAbundances[0] * Zfrac, eps);
        const double g = (mTotSol - mNSol * a - mCSol * b) / mRestSol;
        if (g > eps)
        {
            dN = std::log10(a / g);
            dC = std::log10(b / g);
        }
        else
        {
            // Imported N+C exceed Gutkin's reach (M_N*a + M_C*b >= M_tot). No recipe
            // (dN, dC) can reproduce them. Fall back to the physical log-ratio so the
            // stab lookup stays finite; cell is outside the validated regime.
            dN = std::log10(a);
            dC = std::log10(b);
        }
    }
    else
    {
        // SolarScaled mode: dN, dC come straight from the user-set linear scalings.
        dN = std::log10(solarScalingN());
        dC = std::log10(solarScalingC());
    }
}

////////////////////////////////////////////////////////////////////

void DiffuseIonizedGasMix::bracketDeltaAxis(double value, const vector<double>& axis_values,
                                            const vector<int>& axis_deltaIds, int& deltaIdLo, int& deltaIdHi,
                                            double& w) const
{
    const size_t n = axis_values.size();

    // Endpoint clamp: query at or beyond an endpoint -> degenerate bracket, weight irrelevant.
    if (value <= axis_values.front())
    {
        deltaIdLo = axis_deltaIds.front();
        deltaIdHi = axis_deltaIds.front();
        w = 0.0;
        return;
    }
    if (value >= axis_values.back())
    {
        deltaIdLo = axis_deltaIds.back();
        deltaIdHi = axis_deltaIds.back();
        w = 0.0;
        return;
    }

    // Locate the upper bracket index (smallest hi with axis_values[hi] > value).
    size_t hi = 1;
    while (hi < n && axis_values[hi] <= value) ++hi;
    const size_t lo = hi - 1;
    const double v0 = axis_values[lo];
    const double v1 = axis_values[hi];
    w = (value - v0) / (v1 - v0);
    deltaIdLo = axis_deltaIds[lo];
    deltaIdHi = axis_deltaIds[hi];
}

////////////////////////////////////////////////////////////////////
