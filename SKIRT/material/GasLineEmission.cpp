/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GasLineEmission.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "TextInFile.hpp"
#include <set>

//////////////////////////////////////////////////////////////////////

namespace
{
    // Rest-frame wavelengths [m]
    constexpr double lineWavelengths[GasLineEmission::numLines] = {
        1215.670e-10,  // Lya
        6562.800e-10,  // Ha
        4861.330e-10,  // Hb
        4340.460e-10,  // Hg
        4101.730e-10,  // Hd
        3970.070e-10,  // He
        1.87510e-6,    // Paa
        1.28180e-6,    // Pab
        4.05120e-6,    // Bra
        5875.620e-10,  // He I 5876
        6678.150e-10,  // He I 6678
        7065.190e-10,  // He I 7065
        1.08303e-6,    // He I 10830
        1640.420e-10,  // He II 1640
        4685.680e-10,  // He II 4686
        6548.050e-10,  // [NII] 6548
        6583.460e-10,  // [NII] 6583
        6300.304e-10,  // [OI] 6300
        6363.776e-10,  // [OI] 6364
        3728.815e-10,  // [OII] 3729
        3726.032e-10,  // [OII] 3726
        4363.210e-10,  // [OIII] 4363
        4958.911e-10,  // [OIII] 4959
        5006.843e-10,  // [OIII] 5007
        6716.440e-10,  // [SII] 6716
        6730.810e-10,  // [SII] 6731
    };

    // Particle mass for Doppler broadening [kg]: m_p for H, 4 m_p for He, approximate for metals.
    constexpr double lineMasses[GasLineEmission::numLines] = {
        1.67262192e-27,         // Lya (H)
        1.67262192e-27,         // Ha
        1.67262192e-27,         // Hb
        1.67262192e-27,         // Hg
        1.67262192e-27,         // Hd
        1.67262192e-27,         // He
        1.67262192e-27,         // Paa
        1.67262192e-27,         // Pab
        1.67262192e-27,         // Bra
        4.0 * 1.67262192e-27,   // HeI 5876
        4.0 * 1.67262192e-27,   // HeI 6678
        4.0 * 1.67262192e-27,   // HeI 7065
        4.0 * 1.67262192e-27,   // HeI 10830
        4.0 * 1.67262192e-27,   // HeII 1640
        4.0 * 1.67262192e-27,   // HeII 4686
        14.0 * 1.67262192e-27,  // NII
        14.0 * 1.67262192e-27,  // NII
        16.0 * 1.67262192e-27,  // OI
        16.0 * 1.67262192e-27,  // OI
        16.0 * 1.67262192e-27,  // OII
        16.0 * 1.67262192e-27,  // OII
        16.0 * 1.67262192e-27,  // OIII
        16.0 * 1.67262192e-27,  // OIII
        16.0 * 1.67262192e-27,  // OIII
        32.0 * 1.67262192e-27,  // SII
        32.0 * 1.67262192e-27,  // SII
    };

    // Ion stage index in PhotoIonizationSolver::ionFracs[] of the carrier ion of each collisional line
    // N: stages 12-19 (NI=12, NII=13, ...)
    // O: stages 20-28 (OI=20, OII=21, OIII=22, ...)
    // S: stages 62-66 (SI=62, SII=63, ...)
    constexpr int lineCarrierIonIndex[GasLineEmission::numLines] = {
        -1, -1, -1, -1, -1, -1, -1, -1, -1,  // H lines: not from ion fractions
        -1, -1, -1, -1, -1, -1,              // He lines: handled with the recombination family
        13,                                  // [NII]  -> NII  = ion index 13
        13,                                  // [NII]  -> NII
        20,                                  // [OI]   -> OI   = ion index 20
        20,                                  // [OI]   -> OI
        21,                                  // [OII]  -> OII  = ion index 21
        21,                                  // [OII]  -> OII
        22,                                  // [OIII] -> OIII = ion index 22
        22,                                  // [OIII] -> OIII
        22,                                  // [OIII] -> OIII
        63,                                  // [SII]  -> SII  = ion index 63
        63,                                  // [SII]  -> SII
    };

    // Element index of the carrier of each collisional line (for the abundance lookup)
    // 0=C, 1=N, 2=O, 3=Ne, 4=Mg, 5=Si, 6=S, 7=Fe
    constexpr int lineElementIndex[GasLineEmission::numLines] = {
        -1, -1, -1, -1, -1, -1, -1, -1, -1,  // H lines
        -1, -1, -1, -1, -1, -1,              // He lines
        1,  1,                               // NII -> N
        2,  2,                               // OI -> O
        2,  2,                               // OII -> O
        2,  2,  2,                           // OIII -> O
        6,  6,                               // SII -> S
    };

    // ============== Lyman-alpha recombination data: P_B table (Storey & Hummer 1995) ==============

    // Temperature grid: 10 points from 500 K to 30000 K (non-uniform in log space)
    constexpr int pbNumT = 10;
    constexpr int pbNumNe = 7;

    // log10(T) values computed at compile time from the actual temperatures
    // T = {500, 1000, 3000, 5000, 7500, 10000, 12500, 15000, 20000, 30000}
    // n_e = {1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8}
    inline double logT_SH95(int i)
    {
        static const double v[pbNumT] = {2.69897000, 3.0,        3.47712125, 3.69897000, 3.87506126,
                                         4.0,        4.09691001, 4.17609126, 4.30103000, 4.47712125};
        return v[i];
    }

    // P_B table: probability of emitting a Lyman-alpha photon per recombination
    // Indexed as [temperature][ne], 10x7
    constexpr double P_Lya[pbNumT][pbNumNe] = {
        {0.80793637, 0.90215197, 0.98157401, 0.99705340, 0.99933853, 0.99934845, 0.99967896},
        {0.78439303, 0.87743841, 0.97423464, 0.99680114, 0.99960243, 0.99967816, 0.99976710},
        {0.74159266, 0.82809599, 0.95734434, 0.99359891, 0.99892162, 0.99927343, 0.99979471},
        {0.71946081, 0.80133228, 0.94581237, 0.99246182, 0.99848958, 0.99940329, 0.99949246},
        {0.70117195, 0.77884530, 0.93441344, 0.99075058, 0.99831301, 0.99960223, 0.99940949},
        {0.68808416, 0.76220739, 0.92503048, 0.98912680, 0.99828523, 0.99940740, 0.99944339},
        {0.67755810, 0.74868826, 0.91762730, 0.98784252, 0.99835238, 0.99903921, 0.99977064},
        {0.66862468, 0.73785038, 0.91033402, 0.98701405, 0.99817046, 0.99892625, 0.99953078},
        {0.65534642, 0.72080463, 0.89823705, 0.98497224, 0.99798854, 0.99881330, 0.99929092},
        {0.63566376, 0.69600394, 0.87952243, 0.98121463, 0.99716593, 0.99896144, 0.99905958}};

    /** Bilinear interpolation of P_B(T, n_e) from a Storey & Hummer (1995) table.
        Returns probability of emitting the line photon per Case B recombination.
        T in K, ne in cm^-3. */
    inline double interpolate_PB(const double (&Ptab)[pbNumT][pbNumNe], double T, double ne)
    {
        // Temperature interpolation indices (log space, non-uniform grid)
        int iT = 8;
        double fracT = 1.0;
        if (T <= 500.)
        {
            iT = 0;
            fracT = 0.;
        }
        else if (T < 30000.)
        {
            double logT = std::log10(T);
            iT = 8;
            while (logT < logT_SH95(iT)) iT--;
            fracT = (logT - logT_SH95(iT)) / (logT_SH95(iT + 1) - logT_SH95(iT));
        }

        // Density interpolation indices (log space, uniform 1 dex spacing from 10^2 to 10^8)
        int iN = 0;
        double fracN = 0.;
        if (ne >= 1e8)
        {
            iN = 5;
            fracN = 1.;
        }
        else if (ne > 1e2)
        {
            double d = std::log10(ne) - 2.0;
            double f = std::floor(d);
            fracN = d - f;
            iN = static_cast<int>(f);
        }

        // Bilinear interpolation
        double fLT = 1. - fracT, fLN = 1. - fracN;
        int jT = iT + 1, jN = iN + 1;
        return fLT * (Ptab[iT][iN] * fLN + Ptab[iT][jN] * fracN) + fracT * (Ptab[jT][iN] * fLN + Ptab[jT][jN] * fracN);
    }

    // ============== Physical constants (SI, P_B form) ==============

    // used by the P_B form; Constants::h() differs by 8.8e-8 relative
    constexpr double h_SI = 6.62607015e-34;  // [J s]
    constexpr double c_SI = 2.99792458e8;    // [m/s]

    // ============== Case B recombination-line emissivity tables ==============

    // the built-in recombination lines, mapped to their species cube (0 H I, 1 He I, 2 He II) and
    // their transition. H I and He II are located in the cube by upper/lower quantum number; He I has
    // no quantum numbers in its source table, so it is located by wavelength (upper = lower = -1). No
    // Case B table exists for Lyman-alpha, which uses the P_B form.
    struct LineFile
    {
        int lineIdx;
        int cubeId;
        int upper;
        int lower;
    };
    const vector<LineFile>& recombLineFiles()
    {
        static const vector<LineFile> table = {
            {GasLineEmission::Ha, 0, 3, 2},        {GasLineEmission::Hb, 0, 4, 2},
            {GasLineEmission::Hg, 0, 5, 2},        {GasLineEmission::Hd, 0, 6, 2},
            {GasLineEmission::HeBalmer, 0, 7, 2},  {GasLineEmission::Paa, 0, 4, 3},
            {GasLineEmission::Pab, 0, 5, 3},       {GasLineEmission::Bra, 0, 5, 4},
            {GasLineEmission::HeI5876, 1, -1, -1}, {GasLineEmission::HeI6678, 1, -1, -1},
            {GasLineEmission::HeI7065, 1, -1, -1}, {GasLineEmission::HeI10830, 1, -1, -1},
            {GasLineEmission::HeII1640, 2, 3, 2},  {GasLineEmission::HeII4686, 2, 4, 3}};
        return table;
    }

    // ============== Level populations (statistical equilibrium) ==============

    // clamped log-log interpolation replicating NR::clampedValue<NR::interpolateLogLog>
    double clampedLogLog(double x, const vector<double>& xv, const vector<double>& yv)
    {
        int n = static_cast<int>(xv.size());
        if (x < xv[0]) return yv[0];
        if (x >= xv[n - 1]) return yv[n - 1];
        int i = static_cast<int>(std::upper_bound(xv.begin(), xv.end(), x) - xv.begin()) - 1;
        double x1 = xv[i], x2 = xv[i + 1], f1 = yv[i], f2 = yv[i + 1];
        if (f1 <= 0 || f2 <= 0)
        {
            if (x == x1) return f1;
            if (x == x2) return f2;
            return 0.;
        }
        return f1 * exp(log(x / x1) / log(x2 / x1) * (log(f2 / f1)));
    }

    // LU solve with row pivoting only on a zero diagonal; throws on a singular matrix or a
    // non-finite solution
    vector<double> solveMatrixEquationLU(vector<vector<double>>& matrix)
    {
        size_t size = matrix.size();
        vector<double> solution(size);
        for (size_t i = 0; i < size; i++) solution[i] = matrix[i][size];

        for (size_t k = 0; k < size - 1; ++k)
        {
            if (matrix[k][k] == 0.0)
            {
                size_t maxRow = k;
                double maxAbsVal = std::abs(matrix[k][k]);
                for (size_t i = k + 1; i < size; ++i)
                {
                    double absVal = std::abs(matrix[i][k]);
                    if (absVal > maxAbsVal)
                    {
                        maxAbsVal = absVal;
                        maxRow = i;
                    }
                }
                if (maxRow != k)
                {
                    std::swap(matrix[k], matrix[maxRow]);
                    std::swap(solution[k], solution[maxRow]);
                }
                else
                    throw FATALERROR("Level population matrix is singular at row " + std::to_string(k));
            }
            for (size_t i = k + 1; i < size; ++i)
            {
                double inverse = matrix[i][k] / matrix[k][k];
                for (size_t j = k + 1; j < size; ++j) matrix[i][j] -= inverse * matrix[k][j];
                matrix[i][k] = inverse;
            }
        }
        for (size_t i = 0; i < size; ++i)
            for (size_t j = 0; j < i; ++j) solution[i] -= matrix[i][j] * solution[j];
        for (int i = static_cast<int>(size) - 1; i >= 0; --i)
        {
            for (size_t j = i + 1; j < size; ++j) solution[i] -= matrix[i][j] * solution[j];
            solution[i] /= matrix[i][i];
            if (!std::isfinite(solution[i]))
                throw FATALERROR("Level population solution is not finite at index " + std::to_string(i));
        }
        return solution;
    }

    // one Case B recombination table set for the extended inventory: the species cube index and
    // the line family (used for the carrier mass and the H I / He I / He II flags)
    struct RecombSet
    {
        int cubeId;
        bool isHeI;
        bool isHeII;
        double mass;
    };
    constexpr double _protonMass = 1.67262192e-27;
    const RecombSet _recombSets[] = {
        {0, false, false, _protonMass}, {1, true, false, 4.0 * _protonMass}, {2, false, true, 4.0 * _protonMass}};
}

//////////////////////////////////////////////////////////////////////

GasLineEmission::GasLineEmission() {}

//////////////////////////////////////////////////////////////////////

GasLineEmission::~GasLineEmission() {}

//////////////////////////////////////////////////////////////////////

void GasLineEmission::initialize(const SimulationItem* item)
{
    // remember a pointer to the item and the logger for this simulation
    _item = item;
    _log = item->find<Log>();

    // add the built-in lines to the registry
    _registry.reserve(numLines);
    for (int k = 0; k != numLines; ++k)
    {
        _registry.push_back({lineWavelengths[k], lineMasses[k], lineCarrierIonIndex[k], lineElementIndex[k],
                             k >= HeI5876 && k <= HeI10830, k >= HeII1640 && k <= HeII4686});
        _recombRegistry.table.emplace_back();
        _atomicRegistry.lineModel.push_back(-1);
        _atomicRegistry.lineTransition.push_back(-1);
    }
}

//////////////////////////////////////////////////////////////////////

void GasLineEmission::initializeRecombinationTables()
{
    if (_recombRegistry.ready) return;

    static const char* cubeFiles[3] = {"HI_CaseB_all_lines.stab", "HeI_CaseB_all_lines.stab",
                                       "HeII_CaseB_all_lines.stab"};
    static const char* mapFiles[3] = {"HI_line_index_map.txt", "HeI_line_index_map.txt", "HeII_line_index_map.txt"};

    // open each species' aggregated emissivity cube and read its transition index map
    for (int c = 0; c != 3; ++c)
    {
        _recombRegistry.cube[c].open(_item, cubeFiles[c], "line(1),T(K),n(1/cm3)", "Emis(1)");

        TextInFile in(_item, mapFiles[c], "recombination line index map", true);
        in.addColumn("line index");
        in.addColumn("upper level");
        in.addColumn("lower level");
        in.addColumn("wavelength", "wavelength", "micron");
        double index, upper, lower, wav_m;
        while (in.readRow(index, upper, lower, wav_m))
            _recombRegistry.map[c].push_back(
                {static_cast<int>(index), static_cast<int>(upper), static_cast<int>(lower), wav_m});
    }

    // point each built-in recombination line at its cube line: H I and He II by quantum number,
    // He I by nearest wavelength (its source table carries no quantum numbers)
    for (const auto& entry : recombLineFiles())
    {
        int found = -1;
        if (entry.upper > 0)
        {
            for (const auto& row : _recombRegistry.map[entry.cubeId])
                if (row.upper == entry.upper && row.lower == entry.lower)
                {
                    found = row.index;
                    break;
                }
        }
        else
        {
            double best = std::numeric_limits<double>::max();
            for (const auto& row : _recombRegistry.map[entry.cubeId])
            {
                double dist = std::abs(row.wav_m - lineWavelengths[entry.lineIdx]);
                if (dist < best)
                {
                    best = dist;
                    found = row.index;
                }
            }
            if (best > 5e-10) found = -1;
        }
        if (found < 0)
            throw FATALERROR("No Case B emissivity table for built-in recombination line index "
                             + std::to_string(entry.lineIdx));
        _recombRegistry.table[entry.lineIdx] = {entry.cubeId, found};
    }

    _recombRegistry.ready = true;
}

//////////////////////////////////////////////////////////////////////

void GasLineEmission::initializeAtomicModels()
{
    if (_atomicRegistry.ready) return;

    // carrier species and the lines they serve (LineIndex values)
    struct SpeciesLines
    {
        const char* name;
        vector<int> lines;
    };
    vector<SpeciesLines> table = {{"N_II", {NII6548, NII6583}},
                                  {"O_I", {OI6300, OI6364}},
                                  {"O_II", {OII3729, OII3726}},
                                  {"O_III", {OIII4363, OIII4959, OIII5007}},
                                  {"S_II", {SII6716, SII6731}}};

    for (int k = 0; k != numLines; ++k)
    {
        _atomicRegistry.lineModel[k] = -1;
        _atomicRegistry.lineTransition[k] = -1;
    }
    for (const auto& entry : table)
    {
        AtomicModel model;
        loadAtomicModel(entry.name, {"e-"}, 20, model);
        int slot = static_cast<int>(_atomicRegistry.models.size());
        _atomicRegistry.models.push_back(std::move(model));
        _atomicRegistry.modelNames.push_back(entry.name);
        const auto& m = _atomicRegistry.models[slot];
        for (int lineIdx : entry.lines)
        {
            // map the line to the nearest-wavelength transition of its carrier
            double target = lineWavelengths[lineIdx];
            int kbest = -1;
            double dbest = std::numeric_limits<double>::infinity();
            for (int k = 0; k != m.numLines(); ++k)
            {
                double d = std::abs(m.center[k] - target);
                if (d < dbest)
                {
                    dbest = d;
                    kbest = k;
                }
            }
            if (kbest < 0)
                throw FATALERROR("Atomic model " + string(entry.name)
                                 + " has no radiative transition for built-in collisional line index "
                                 + std::to_string(lineIdx));
            _atomicRegistry.lineModel[lineIdx] = slot;
            _atomicRegistry.lineTransition[lineIdx] = kbest;
        }
    }
    _atomicRegistry.ready = true;
}

//////////////////////////////////////////////////////////////////////

void GasLineEmission::initializeExtendedLineRegistry(const vector<SpeciesSpec>& species)
{
    // ensure that the recombination tables and atomic models are loaded
    initializeRecombinationTables();
    initializeAtomicModels();

    // quit if extended lines are already loaded
    if (static_cast<int>(_registry.size()) > numLines) return;

    // (cube, line index) pairs already claimed by the built-in recombination lines
    std::set<std::pair<int, int>> claimedCubeLines;
    for (const auto& entry : recombLineFiles())
        claimedCubeLines.insert(
            {_recombRegistry.table[entry.lineIdx].cubeId, _recombRegistry.table[entry.lineIdx].lineIdx});

    // (model slot, transition) pairs already claimed by the built-in collisional lines
    std::set<std::pair<int, int>> claimedTransitions;
    for (int k = 0; k != numLines; ++k)
        if (_atomicRegistry.lineModel[k] >= 0)
            claimedTransitions.insert({_atomicRegistry.lineModel[k], _atomicRegistry.lineTransition[k]});

    int numAdded = 0;
    int numModelsLoaded = 0;

    // recombination inventory: every line in each species' cube, enumerated from its index map
    for (const auto& set : _recombSets)
    {
        for (const auto& row : _recombRegistry.map[set.cubeId])
        {
            if (!claimedCubeLines.insert({set.cubeId, row.index}).second) continue;

            _recombRegistry.table.push_back({set.cubeId, row.index});
            _registry.push_back({row.wav_m, set.mass, -1, -1, set.isHeI, set.isHeII});
            _atomicRegistry.lineModel.push_back(-1);
            _atomicRegistry.lineTransition.push_back(-1);
            ++numAdded;
        }
    }

    // collisional lines: every transition of every loadable species
    for (const auto& spec : species)
    {
        int slot = -1;
        for (int m = 0; m != static_cast<int>(_atomicRegistry.modelNames.size()); ++m)
            if (_atomicRegistry.modelNames[m] == spec.name) slot = m;
        if (slot < 0)
        {
            // not every candidate (element, ionization-stage) species has data in the resource
            // pack; skip absent ones by checking first rather than treating that as an error
            if (!FilePaths::hasResource(spec.name + "_Mass.txt")) continue;

            AtomicModel model;
            loadAtomicModel(spec.name, {"e-"}, 20, model);
            ++numModelsLoaded;
            slot = static_cast<int>(_atomicRegistry.models.size());
            _atomicRegistry.models.push_back(std::move(model));
            _atomicRegistry.modelNames.push_back(spec.name);
        }
        const auto& m = _atomicRegistry.models[slot];
        for (int t = 0; t != m.numLines(); ++t)
        {
            if (claimedTransitions.count({slot, t})) continue;
            _registry.push_back({m.center[t], m.mass, spec.carrierIonIndex, spec.elementIndex, false, false});
            _recombRegistry.table.emplace_back();
            _atomicRegistry.lineModel.push_back(slot);
            _atomicRegistry.lineTransition.push_back(t);
            ++numAdded;
        }
    }

    // inform the user
    _log->info("Line registry extended with " + std::to_string(numAdded) + " lines from "
               + std::to_string(numModelsLoaded) + " atomic models");
}

//////////////////////////////////////////////////////////////////////

const vector<GasLineEmission::LineDef>& GasLineEmission::lineRegistry() const
{
    return _registry;
}

//////////////////////////////////////////////////////////////////////

double GasLineEmission::hydrogenLineLuminosity(double T, double ne, double gammaHI, double nHI, double V_cm3) const
{
    double PB = interpolate_PB(P_Lya, T, ne);
    double lambda_m = lineWavelengths[Lya];
    double hnu = h_SI * c_SI / lambda_m;  // photon energy [J]
    // L = h*nu * P_B * Gamma_HI * n_HI * V [W]
    return hnu * PB * gammaHI * nHI * V_cm3;
}

//////////////////////////////////////////////////////////////////////

double GasLineEmission::recombinationLineLuminosity(int lineIdx, double T, double ne, double nIon, double gammaHI,
                                                    double nHI, double V_cm3) const
{
    // Lyman-alpha has no Case B table and uses the P_B form
    if (lineIdx == Lya) return hydrogenLineLuminosity(T, ne, gammaHI, nHI, V_cm3);

    const auto& ref = _recombRegistry.table[lineIdx];
    double eps = _recombRegistry.cube[ref.cubeId](static_cast<double>(ref.lineIdx), T, ne);
    // eps [W m3] with cgs densities and volume: L [W] = eps * ne * nIon * V_cm3 * 1e6
    return eps * ne * nIon * V_cm3 * 1e6;
}

//////////////////////////////////////////////////////////////////////

int GasLineEmission::lineModelSlot(int lineIdx) const
{
    if (lineIdx < 0 || lineIdx >= static_cast<int>(_atomicRegistry.lineModel.size())) return -1;
    return _atomicRegistry.lineModel[lineIdx];
}

//////////////////////////////////////////////////////////////////////

int GasLineEmission::lineTransition(int lineIdx) const
{
    if (lineIdx < 0 || lineIdx >= static_cast<int>(_atomicRegistry.lineTransition.size())) return -1;
    return _atomicRegistry.lineTransition[lineIdx];
}

//////////////////////////////////////////////////////////////////////

void GasLineEmission::loadAtomicModel(const string& speciesName, const vector<string>& partnerNames, int maxNumLevels,
                                      AtomicModel& model) const
{
    model = AtomicModel();

    // mass
    {
        TextInFile infile(_item, speciesName + "_Mass.txt", "mass", true, true);
        infile.addColumn("Mass", "mass", "amu");
        double mass;
        if (infile.readRow(mass)) model.mass = mass;
    }

    // energy levels and statistical weights, capped at maxNumLevels
    {
        TextInFile infile(_item, speciesName + "_Energy.txt", "energy levels", true, true);
        infile.addColumn("Energy", "energy", "1/cm");
        infile.addColumn("Weight");
        double energy, weight;
        while (infile.readRow(energy, weight))
        {
            model.energy.push_back(energy);
            model.weight.push_back(weight);
            if (static_cast<int>(model.energy.size()) == maxNumLevels) break;
        }
    }
    int numLevels = model.numLevels();

    // radiative transitions, dropping any that involve levels beyond the cap
    {
        TextInFile infile(_item, speciesName + "_Rad_Coeff.txt", "radiative transitions", true, true);
        infile.addColumn("Up index");
        infile.addColumn("Low index");
        infile.addColumn("Einstein A", "transitionrate", "1/s");
        double up, low, rate;
        while (infile.readRow(up, low, rate))
        {
            int indexUp = static_cast<int>(up), indexLow = static_cast<int>(low);
            if (indexUp < numLevels && indexLow < numLevels)
            {
                model.indexUpRad.push_back(indexUp);
                model.indexLowRad.push_back(indexLow);
                model.einsteinA.push_back(rate);
            }
        }
    }
    int numRad = model.numLines();

    // derived per-line data: centers, Einstein B coefficients, branching ratios
    model.center.resize(numRad);
    model.einsteinBul.resize(numRad);
    model.einsteinBlu.resize(numRad);
    model.branchRatio.resize(numRad);
    vector<double> sumA(numLevels, 0.);
    for (int k = 0; k != numRad; ++k) sumA[model.indexUpRad[k]] += model.einsteinA[k];
    for (int k = 0; k != numRad; ++k)
    {
        int up = model.indexUpRad[k], low = model.indexLowRad[k];
        double lambda = Constants::h() * Constants::c() / (model.energy[up] - model.energy[low]);
        model.center[k] = lambda;
        model.einsteinBul[k] =
            model.einsteinA[k] * pow(lambda, 5.) / (2. * Constants::h() * Constants::c() * Constants::c());
        model.einsteinBlu[k] = model.einsteinBul[k] * model.weight[up] / model.weight[low];
        model.branchRatio[k] = sumA[up] > 0. ? model.einsteinA[k] / sumA[up] : 0.;
    }

    // collisional transitions per partner
    for (const auto& pname : partnerNames)
    {
        ColPartner partner;
        partner.name = pname;

        {
            TextInFile infile(_item, speciesName + "_Col_" + pname + "_Temp.txt", "temperature grid", true, true);
            infile.addColumn("Temperature", "temperature", "K");
            Array tgrid;
            infile.readAllColumns(tgrid);
            partner.T.assign(begin(tgrid), end(tgrid));
        }
        {
            int numTemperatures = static_cast<int>(partner.T.size());
            TextInFile infile(_item, speciesName + "_Col_" + pname + "_Coeff.txt", "collisional transitions", true,
                              true);
            infile.addColumn("Up index");
            infile.addColumn("Low index");
            for (int i = 0; i != numTemperatures; ++i) infile.addColumn("Collisional K", "collisionalrate", "cm3/s");
            Array row;
            while (infile.readRow(row))
            {
                int up = static_cast<int>(row[0]), low = static_cast<int>(row[1]);
                if (up < numLevels && low < numLevels)
                {
                    partner.indexUpCol.push_back(up);
                    partner.indexLowCol.push_back(low);
                    vector<double> K(numTemperatures);
                    for (int i = 0; i != numTemperatures; ++i) K[i] = row[i + 2];
                    partner.Kul.push_back(std::move(K));
                }
            }
        }
        model.colPartner.push_back(std::move(partner));
    }
}

//////////////////////////////////////////////////////////////////////

const GasLineEmission::AtomicModel& GasLineEmission::atomicModel(int slot) const
{
    return _atomicRegistry.models[slot];
}

//////////////////////////////////////////////////////////////////////

vector<double> GasLineEmission::solveLevelPopulations(const AtomicModel& model, const Environment& env) const
{
    int numLevels = model.numLevels();
    vector<vector<double>> matrix(numLevels, vector<double>(numLevels + 1, 0.));

    // radiative transitions: spontaneous emission always; pumping terms only if meanJ is provided
    bool pumping = !env.meanJ.empty();
    for (int k = 0; k != model.numLines(); ++k)
    {
        int up = model.indexUpRad[k], low = model.indexLowRad[k];
        matrix[up][up] -= model.einsteinA[k];
        matrix[low][up] += model.einsteinA[k];
        if (pumping)
        {
            double J = env.meanJ[k];
            matrix[up][up] -= model.einsteinBul[k] * J;
            matrix[low][up] += model.einsteinBul[k] * J;
            matrix[low][low] -= model.einsteinBlu[k] * J;
            matrix[up][low] += model.einsteinBlu[k] * J;
        }
    }

    // collisional transitions; Kul, Klu and the partner density are floored to keep the matrix regular
    double T = env.Tkin;
    for (size_t c = 0; c != model.colPartner.size(); ++c)
    {
        const auto& partner = model.colPartner[c];
        double Tmin = partner.T.front(), Tmax = partner.T.back();
        for (size_t t = 0; t != partner.indexUpCol.size(); ++t)
        {
            int up = partner.indexUpCol[t], low = partner.indexLowCol[t];
            double weightRatio = model.weight[up] / model.weight[low];
            double energyDiff = model.energy[up] - model.energy[low];
            double Kconversion = std::max(weightRatio * exp(-energyDiff / Constants::k() / T), 1e-15);
            double Trep = std::max(Tmin, std::min(T, Tmax));
            double Kul = clampedLogLog(Trep, partner.T, partner.Kul[t]);
            double Klu = Kul * Kconversion;
            if (Kul <= 0.)
            {
                _log->warning("collisional rate Kul <= 0 for partner " + partner.name + " transition ("
                              + std::to_string(up) + "-" + std::to_string(low) + "); flooring to 1e-20");
                Kul = 1.0e-20;
            }
            if (Klu <= 0.)
            {
                _log->warning("collisional rate Klu <= 0 for partner " + partner.name + " transition ("
                              + std::to_string(up) + "-" + std::to_string(low) + "); flooring to 1e-20*Kul");
                Klu = 1.0e-20 * Kul;
            }
            double n = std::max(env.nPartner[c], 1.0e-20);
            matrix[up][up] -= Kul * n;
            matrix[low][low] -= Klu * n;
            matrix[up][low] += Klu * n;
            matrix[low][up] += Kul * n;
        }
    }

    // replace the last row by the number density normalization and solve
    for (int p = 0; p != numLevels; ++p) matrix[numLevels - 1][p] = 1.;
    matrix[numLevels - 1][numLevels] = env.nTotal;
    return solveMatrixEquationLU(matrix);
}

//////////////////////////////////////////////////////////////////////

vector<double> GasLineEmission::lineEmissivities(const AtomicModel& model, const vector<double>& pops) const
{
    vector<double> eps(model.numLines());
    for (int k = 0; k != model.numLines(); ++k)
        eps[k] = Constants::h() * Constants::c() / model.center[k] * model.einsteinA[k] * pops[model.indexUpRad[k]];
    return eps;
}

//////////////////////////////////////////////////////////////////////
