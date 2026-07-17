/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PhotoIonizationSolver.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"
#include "PhotoIonizationRates.hpp"
#include "VernerCrossSections.hpp"
#include <cstring>

//////////////////////////////////////////////////////////////////////

// Out-of-line definitions for static constexpr members (required for ODR-use on GCC/Linux)
constexpr int PhotoIonizationSolver::numStages[];
constexpr int PhotoIonizationSolver::stageOffset[];

//////////////////////////////////////////////////////////////////////

namespace
{
    // physical constants in CGS
    constexpr double eV_erg = 1.60217663e-12;  // 1 eV in erg
    constexpr double kB_cgs = 1.380649e-16;    // Boltzmann constant [erg/K]
    constexpr double fourPi = 4. * M_PI;

    // convergence parameters for electron density iteration
    constexpr int maxNeIterations = 100;
    constexpr double neTolerance = 1e-9;  // absolute tolerance on xe
    constexpr double minElectronFraction = 1e-10;

    // cosmic-ray background heating constants
    // Indriolo et al. (2007, ApJ 671, 1736): zeta_CR ~ 2e-16 s^-1 per H atom
    // Wolfire et al. (1995, ApJ 443, 152): ~35 eV deposited as heat per CR ionization
    constexpr double zetaCR = 2.0e-16;                             // CR ionization rate [s^-1]
    constexpr double crHeatPerIonization = 35.0 * 1.60217663e-12;  // 35 eV [erg]

    // temperature bisection parameters
    constexpr double Tmin = 300.;
    constexpr double Tmax = 1e9;
    constexpr int maxTIterations = 100;
    constexpr double Ttolerance = 1e-3;  // relative tolerance on T

    // number of ionization stages per element (including fully ionized)
    // H: HI, HII (2); He: HeI, HeII, HeIII (3); C: CI-CVII (7); N: NI-NVIII (8);
    // O: OI-OIX (9); Ne: NeI-NeIX (9); Mg: MgI-MgXI (11); Si: SiI-SiXIII (13);
    // S: SI-SV (5); Fe: FeI-FeXVII (17)
    constexpr int nStages[] = {2, 3, 7, 8, 9, 9, 11, 13, 5, 17};
    constexpr int nElements = 10;

    // offset of first ionization stage in the ion fraction array, for each element
    constexpr int offset[] = {0, 2, 5, 12, 20, 29, 38, 49, 62, 67};

    // offset of the first ion of each element in the Verner cross-section indexing
    // H=0, He=1, C=3, N=9, O=16, Ne=24, Mg=32, Si=42, S=54, Fe=58
    constexpr int vernerOffset[] = {0, 1, 3, 9, 16, 24, 32, 42, 54, 58};

    // number of ions with cross-sections per element (excludes fully ionized)
    constexpr int nIonsXS[] = {1, 2, 6, 7, 8, 8, 10, 12, 4, 16};

    // safe division: a/b, returns 0 if b ~ 0
    inline double safediv(double a, double b)
    {
        return (b > 1e-300) ? a / b : 0.;
    }

    // Total thermally-averaged free-free Gaunt factor from van Hoof et al. (2014)
    // gamma2 = Z^2 * Ry / (kT) where Ry = 2.1799e-11 erg = 13.606 eV
    inline double totalGauntFactor(double gamma2)
    {
        if (gamma2 <= 1e-6)
        {
            double g = sqrt(gamma2);
            return 1.102635 + 1.186 * g + 0.86 * g * g;
        }
        if (gamma2 >= 1e10)
        {
            double g = sqrt(gamma2);
            return 1. + pow(g, -2. / 3.);
        }
        double logG = log10(gamma2);
        if (logG <= 0.8)
        {
            double l2 = logG * logG, l3 = l2 * logG, l4 = l3 * logG;
            double num = 1.43251926625281 + 3.50626935257777e-1 * logG + 4.36183448595035e-1 * l2
                         + 6.03536387105599e-2 * l3 + 3.66626405363100e-2 * l4;
            double den = 1. + 2.92525161994346e-1 * logG + 4.05566949766954e-1 * l2 + 5.62573012783879e-2 * l3
                         + 3.33019373823972e-2 * l4;
            return num / den;
        }
        else
        {
            double l2 = logG * logG, l3 = l2 * logG, l4 = l3 * logG;
            double num = 1.45481634667278 - 9.55399384620923e-2 * logG + 1.46327814151538e-1 * l2
                         - 1.41489406498468e-2 * l3 + 2.74891413242655e-3 * l4;
            double den = 1. + 3.31149751183539e-2 * logG + 1.31127367293310e-1 * l2 - 1.32658217746618e-2 * l3
                         + 2.74809263365693e-3 * l4;
            return num / den;
        }
    }
}

//////////////////////////////////////////////////////////////////////

void PhotoIonizationSolver::initialize(const Array& lambdav, const Array& dlambdav, bool coolingTable)
{
    _numBins = lambdav.size();
    _lambda.resize(_numBins);
    _dlambda.resize(_numBins);
    _energy.resize(_numBins);
    _energyErg.resize(_numBins);

    for (int i = 0; i < _numBins; i++)
    {
        _lambda[i] = lambdav[i];
        _dlambda[i] = dlambdav[i];
        // photon energy: E = hc/lambda, convert to eV and erg
        double E_J = Constants::h() * Constants::c() / lambdav[i];  // in Joules
        _energyErg[i] = E_J * 1e7;                                  // J -> erg
        _energy[i] = E_J / Constants::Qelectron();                  // J -> eV
    }

    // load metal cooling table if path provided
    if (coolingTable) loadCoolingTable();

    // default: effective stages = full stages
    for (int e = 0; e < nElements; e++) _nStagesEff[e] = nStages[e];
}

//////////////////////////////////////////////////////////////////////

int PhotoIonizationSolver::totalEffectiveStages() const
{
    int total = 0;
    for (int e = 0; e < nElements; e++) total += _nStagesEff[e];
    return total;
}

//////////////////////////////////////////////////////////////////////

void PhotoIonizationSolver::setMaxIonizationEnergy(double maxEV)
{
    for (int e = 0; e < nElements; e++)
    {
        if (maxEV < 0.)
        {
            _nStagesEff[e] = nStages[e];
            continue;
        }
        int voff = vernerOffset[e];
        int nxs = nIonsXS[e];
        int lastAllowed = -1;
        for (int s = 0; s < nxs; s++)
        {
            if (VernerCrossSections::ionizationPotential(voff + s) <= maxEV) lastAllowed = s;
        }
        // stages 0..lastAllowed + the sink at lastAllowed+1; ensure at least 2
        _nStagesEff[e] = std::max(2, std::min(lastAllowed + 2, nStages[e]));
    }
}

//////////////////////////////////////////////////////////////////////

void PhotoIonizationSolver::setCosmicRayHeating(bool enabled, double scaleFactor)
{
    _cosmicRayHeating = enabled;
    _cosmicRayScale = scaleFactor;
}

//////////////////////////////////////////////////////////////////////

void PhotoIonizationSolver::initializeIonFractions(double ionFracs[]) const
{
    // zero all stages first (ensures unused capped stages are clean)
    for (int i = 0; i < totalStages; i++) ionFracs[i] = 0.;

    // H: mostly ionized
    ionFracs[0] = 1e-4;       // HI
    ionFracs[1] = 1. - 1e-4;  // HII

    // He: singly ionized
    ionFracs[2] = 1e-4;              // HeI
    ionFracs[3] = 1. - 1e-4 - 1e-6;  // HeII
    ionFracs[4] = 1e-6;              // HeIII

    // metals: singly ionized
    for (int e = 2; e < nElements; e++) ionFracs[offset[e] + 1] = 1.;
}

//////////////////////////////////////////////////////////////////////

void PhotoIonizationSolver::computePhotoionizationRates(const Array& Jv, double gamma[]) const
{
    // initialize to zero
    for (int ion = 0; ion < VernerCrossSections::numIons; ion++) gamma[ion] = 0.;

    // integrate 4*pi * J_lambda * sigma(E) / (h*nu) * dlambda for each ion
    // J_lambda [W m^-3 sr^-1], convert to J_nu [erg s^-1 cm^-2 Hz^-1 sr^-1]:
    //   J_nu = J_lambda * lambda^2 / c  [W m^-2 Hz^-1 sr^-1]
    //   -> convert to CGS: * 1e3 [erg s^-1 cm^-2 Hz^-1 sr^-1]
    // photoionization rate = 4pi * integral{ J_nu * sigma / (h*nu) d_nu }
    //                      = 4pi * integral{ J_lambda * sigma / (h*c/lambda) * dlambda }
    //                      = 4pi * integral{ J_lambda * sigma * lambda / (h*c) * dlambda }
    // With J_lambda in [W m^-3 sr^-1] and sigma in [cm^2]:
    //   gamma [s^-1] = 4pi * sum_i { J_lambda_i * sigma(E_i) * lambda_i / (h*c) * dlambda_i }
    // where h*c in SI = 1.98645e-25 J*m, and J_lambda is W/m^3/sr = J/(s*m^3*sr)
    // So gamma = 4pi * sum { J_lambda * dlambda * lambda * sigma / (hc) }
    // units: [W/m^3/sr] * [m] * [m] * [cm^2] / [J*m] = [1/s] * [cm^2/m^2]
    // need to convert cm^2 to m^2: factor 1e-4, or equivalently divide by 1e4

    double hc_SI = Constants::h() * Constants::c();  // h*c in J*m

    for (int i = 0; i < _numBins; i++)
    {
        double E_eV = _energy[i];
        double lam = _lambda[i];
        double dl = _dlambda[i];
        double Jlam = Jv[i];  // W m^-3 sr^-1
        if (Jlam <= 0.) continue;

        // common factor: 4pi * J_lambda * dlambda * lambda / (hc) [s^-1 per cm^2 -> need * 1e-4 for m^2]
        double factor = fourPi * Jlam * dl * lam / hc_SI * 1e-4;

        for (int ion = 0; ion < VernerCrossSections::numIons; ion++)
        {
            double sig = VernerCrossSections::sigma(ion, E_eV);
            if (sig > 0.) gamma[ion] += factor * sig;
        }
    }
}

//////////////////////////////////////////////////////////////////////

double PhotoIonizationSolver::computeHeatingRate(const Array& Jv, const double ionFracs[], double nH, double yHe,
                                                 const double metalAbundances[8], double ne, const int effNS[]) const
{
    (void)ne;  // reserved for future Compton heating term

    // heating = sum over ions of n_ion * integral{ 4pi * J_nu * sigma(nu) * (h*nu - IP) / (h*nu) dnu }
    // = sum over ions of n_ion * integral{ 4pi * J_lambda * sigma(E) * (E - IP) / E * dlambda * lambda / (hc) }
    // same unit conversion as photoionization rates, but weighted by (E - IP) in erg

    double hc_SI = Constants::h() * Constants::c();
    double heating = 0.;

    // precompute number densities for each ion with XS
    double nIon[VernerCrossSections::numIons];
    for (int i = 0; i < VernerCrossSections::numIons; i++) nIon[i] = 0.;
    // H
    nIon[0] = nH * ionFracs[0];  // n_HI
    // He
    double nHe = nH * yHe;
    nIon[1] = nHe * ionFracs[2];  // n_HeI
    nIon[2] = nHe * ionFracs[3];  // n_HeII
    // metals: element e (starting from C at element index 2), only up to effective stages
    for (int e = 2; e < nElements; e++)
    {
        double nElem = nH * metalAbundances[e - 2];
        int off = offset[e];
        int voff = vernerOffset[e];
        int nxs = effNS[e] - 1;  // effective number of ions with XS
        for (int s = 0; s < nxs; s++) nIon[voff + s] = nElem * ionFracs[off + s];
    }

    for (int i = 0; i < _numBins; i++)
    {
        double E_eV = _energy[i];
        double lam = _lambda[i];
        double dl = _dlambda[i];
        double Jlam = Jv[i];
        if (Jlam <= 0.) continue;

        double baseFactor = fourPi * Jlam * dl * lam / hc_SI * 1e-4;

        for (int ion = 0; ion < VernerCrossSections::numIons; ion++)
        {
            double sig = VernerCrossSections::sigma(ion, E_eV);
            if (sig <= 0.) continue;
            double IP = VernerCrossSections::ionizationPotential(ion);
            if (E_eV <= IP) continue;
            double excessErg = (E_eV - IP) * eV_erg;
            heating += nIon[ion] * baseFactor * sig * excessErg;
        }
    }

    // Cosmic-ray background heating: Gamma_CR = scale * zeta_CR * 35 eV * nHI
    if (_cosmicRayHeating)
    {
        double nHI = nH * ionFracs[0];
        heating += _cosmicRayScale * zetaCR * crHeatPerIonization * nHI;
    }

    return heating;  // [erg s^-1 cm^-3]
}

//////////////////////////////////////////////////////////////////////

namespace
{
    // Parse ion name like "OIII" -> (element, stage) -> ionFracs index.
    // Returns -1 if the ion is out of range of the tracked stages.
    int parseIonName(const std::string& name)
    {
        // Element prefixes ordered longest-first to avoid "N" matching "Ne"
        struct ElemDef
        {
            const char* sym;
            int offset;
            int nstages;
        };
        static const ElemDef elems[] = {{"Ne", 29, 9}, {"Mg", 38, 11}, {"Si", 49, 13}, {"Fe", 67, 17},
                                        {"C", 5, 7},   {"N", 12, 8},   {"O", 20, 9},   {"S", 62, 5}};
        // Roman numerals
        struct RomanDef
        {
            const char* str;
            int val;
        };
        static const RomanDef romans[] = {
            {"XXVI", 26}, {"XXV", 25},   {"XXIV", 24}, {"XXIII", 23}, {"XXII", 22}, {"XXI", 21}, {"XX", 20},
            {"XIX", 19},  {"XVIII", 18}, {"XVII", 17}, {"XVI", 16},   {"XV", 15},   {"XIV", 14}, {"XIII", 13},
            {"XII", 12},  {"XI", 11},    {"X", 10},    {"IX", 9},     {"VIII", 8},  {"VII", 7},  {"VI", 6},
            {"V", 5},     {"IV", 4},     {"III", 3},   {"II", 2},     {"I", 1}};

        for (auto& e : elems)
        {
            size_t slen = std::strlen(e.sym);
            if (name.substr(0, slen) == e.sym)
            {
                std::string roman = name.substr(slen);
                for (auto& r : romans)
                {
                    if (roman == r.str)
                    {
                        int stage = r.val;  // I = neutral = stage 1
                        int idx = e.offset + (stage - 1);
                        if (stage - 1 >= e.nstages) return -1;
                        return idx;
                    }
                }
            }
        }
        return -1;
    }
}

//////////////////////////////////////////////////////////////////////

void PhotoIonizationSolver::loadCoolingTable()
{
    throw FATALERROR("Loading the cooling table is not yet implemented");

    parseIonName("");
    _hasCoolingTable = true;
}

//////////////////////////////////////////////////////////////////////

double PhotoIonizationSolver::interpolateMetalCooling(int ionFracsIdx, double T, double ne) const
{
    if (!_hasCoolingTable) return 0.;
    if (ionFracsIdx < _metalOffset || ionFracsIdx >= totalStages) return 0.;

    int ionSlot = ionFracsIdx - _metalOffset;
    double logT = std::log10(std::max(T, 1.));
    double logNe = std::log10(std::max(ne, 1e-30));

    // Clamp to grid boundaries
    logT = std::max(logT, _coolLogT.front());
    logT = std::min(logT, _coolLogT.back());

    // Below the ne grid minimum, cooling per ion scales linearly with ne
    // (all lines are sub-critical, coronal approximation is exact)
    double neScaleFactor = 1.;
    if (logNe < _coolLogNe.front())
    {
        neScaleFactor = ne / std::pow(10., _coolLogNe.front());
        logNe = _coolLogNe.front();
    }
    else
    {
        logNe = std::min(logNe, _coolLogNe.back());
    }

    // Find bracketing indices for T
    int iT = 0;
    for (int i = 0; i < _coolNT - 1; i++)
    {
        if (_coolLogT[i + 1] >= logT)
        {
            iT = i;
            break;
        }
    }
    if (logT >= _coolLogT.back()) iT = _coolNT - 2;

    // Find bracketing indices for ne
    int iNe = 0;
    for (int i = 0; i < _coolNNe - 1; i++)
    {
        if (_coolLogNe[i + 1] >= logNe)
        {
            iNe = i;
            break;
        }
    }
    if (logNe >= _coolLogNe.back()) iNe = _coolNNe - 2;

    // Fractional positions
    double dT = _coolLogT[iT + 1] - _coolLogT[iT];
    double tT = (dT > 0.) ? (logT - _coolLogT[iT]) / dT : 0.;

    double dNe = _coolLogNe[iNe + 1] - _coolLogNe[iNe];
    double tNe = (dNe > 0.) ? (logNe - _coolLogNe[iNe]) / dNe : 0.;

    // Bilinear interpolation in log space on the cooling values (which can be zero)
    // Use log(Lambda) for interpolation where Lambda > 0, linear otherwise
    auto getVal = [&](int jT, int jNe) -> double {
        return _coolData[ionSlot * (_coolNT * _coolNNe) + jT * _coolNNe + jNe];
    };

    double v00 = getVal(iT, iNe);
    double v10 = getVal(iT + 1, iNe);
    double v01 = getVal(iT, iNe + 1);
    double v11 = getVal(iT + 1, iNe + 1);

    // If all positive, interpolate in log space for better accuracy
    if (v00 > 0. && v10 > 0. && v01 > 0. && v11 > 0.)
    {
        double lv00 = std::log(v00), lv10 = std::log(v10);
        double lv01 = std::log(v01), lv11 = std::log(v11);
        double lv0 = lv00 + tT * (lv10 - lv00);
        double lv1 = lv01 + tT * (lv11 - lv01);
        return neScaleFactor * std::exp(lv0 + tNe * (lv1 - lv0));
    }
    else
    {
        // Linear interpolation (handles zeros)
        double v0 = v00 + tT * (v10 - v00);
        double v1 = v01 + tT * (v11 - v01);
        return neScaleFactor * std::max(v0 + tNe * (v1 - v0), 0.);
    }
}

//////////////////////////////////////////////////////////////////////

double PhotoIonizationSolver::computeCoolingRate(double T, const double ionFracs[], double nH, double yHe,
                                                 const double metalAbundances[8], double ne, const int effNS[]) const
{
    double cooling = 0.;

    // ion number densities
    double nHI = nH * ionFracs[0];
    double nHII = nH * ionFracs[1];
    double nHe = nH * yHe;
    double nHeI = nHe * ionFracs[2];
    double nHeII = nHe * ionFracs[3];
    double nHeIII = nHe * ionFracs[4];

    // -- Hydrogen cooling --

    // Recombination cooling (Case B): Hui & Gnedin (1997) fitted form
    {
        double lambda = 315614. / T;  // 2 * T_HI / T
        cooling += T * 3.435e-30 * pow(lambda, 1.970) / pow(1. + pow(lambda / 2.250, 0.376), 3.720) * ne * nHII;
    }

    // Collisional ionization cooling: beta * IP * ne * nHI
    cooling += PhotoIonizationRates::betaHI(T) * VernerCrossSections::HI_eV * eV_erg * ne * nHI;

    // HI collisional excitation cooling (Lyman series: Ly-alpha through Ly-delta, N=2-5)
    // Maxwellian-averaged effective collision strengths (Anderson et al. 2000, 2002)
    {
        // Transition energies [erg]: E_1n = (1 - 1/n^2) * E_H, where E_H = 13.606 eV
        constexpr double E_H_erg = 13.6057 * 1.60217663e-12;  // hydrogen ionization energy in erg
        constexpr double energyLya = 0.75 * E_H_erg;          // Ly-alpha
        constexpr double energyLyb = (8. / 9.) * E_H_erg;     // Ly-beta
        constexpr double energyLyg = (15. / 16.) * E_H_erg;   // Ly-gamma
        constexpr double energyLyd = (24. / 25.) * E_H_erg;   // Ly-delta
        constexpr double tempLya = energyLya / kB_cgs;        // 118348 K
        constexpr double tempLyb = energyLyb / kB_cgs;        // 140573 K
        constexpr double tempLyg = energyLyg / kB_cgs;        // 148304 K
        constexpr double tempLyd = energyLyd / kB_cgs;        // 151887 K

        // q_coef = h/(4*pi*me) * sqrt(h^3 / (2*pi*me*kB)) ~ 4.3e-6 in CGS
        constexpr double q_coef = 4.339e-6;  // pre-computed coefficient [cm^3 erg / s / K^{1/2}]

        double invSqrtT = 1. / sqrt(T);
        double T6 = 1e-6 * T;
        double HI_cool;
        if (T >= 3e5)
        {
            // High-T regime: constant collision strengths
            HI_cool = (6.107169036302359e-11 * exp(-tempLya / T) + 1.5693136239683913e-11 * exp(-tempLyb / T)
                       + 6.665146531965045e-12 * exp(-tempLyg / T) + 3.4378000369213142e-12 * exp(-tempLyd / T))
                      * q_coef * invSqrtT;
        }
        else
        {
            // Low-T regime: polynomial-fitted collision strengths from Anderson et al.
            HI_cool = (energyLya * (0.616414 + T6 * (16.8152 + T6 * (-32.0571 + 35.5428 * T6))) * exp(-tempLya / T)
                       + energyLyb * (0.217382 + T6 * (3.92604 + T6 * (-10.6349 + 13.7721 * T6))) * exp(-tempLyb / T)
                       + energyLyg * (0.0959324 + T6 * (1.89951 + T6 * (-6.96467 + 10.6362 * T6))) * exp(-tempLyg / T)
                       + energyLyd * (0.0747075 + T6 * (0.670939 + T6 * (-2.28512 + 3.4796 * T6))) * exp(-tempLyd / T))
                      * q_coef * invSqrtT;
        }
        cooling += HI_cool * ne * nHI;
    }

    // -- Helium cooling --

    // HeII recombination cooling
    cooling += PhotoIonizationRates::alphaBHeII(T) * kB_cgs * T * ne * nHeII;

    // HeII dielectronic recombination cooling (0.75 * alpha_DR * T_HeI * kB)
    cooling += 0.75 * PhotoIonizationRates::alphaDRHeII(T) * 631515. * kB_cgs * ne * nHeII;

    // HeIII recombination cooling: Hui & Gnedin (1997) fitted form (Z^2=4 scaling)
    {
        double lambda = 2. * 631515. / T;
        cooling += T * 8. * 3.435e-30 * pow(lambda, 1.970) / pow(1. + pow(lambda / 2.250, 0.376), 3.720) * ne * nHeIII;
    }

    // HeI collisional ionization cooling
    cooling += PhotoIonizationRates::betaHeI(T) * VernerCrossSections::HeI_eV * eV_erg * ne * nHeI;

    // HeII collisional ionization cooling
    cooling += PhotoIonizationRates::betaHeII(T) * VernerCrossSections::HeII_eV * eV_erg * ne * nHeII;

    // HeI collisional excitation cooling (He 2^3S, 19.8 eV)
    cooling += 9.1e-27 * pow(T, -0.1687) * exp(-213751. / T) / (1. + sqrt(T * 1e-5)) * ne * nHeI;

    // HeII collisional excitation cooling (He+ Ly-alpha)
    cooling += 5.54e-17 * pow(T, -0.397) * exp(-473638. / T) / (1. + sqrt(T * 1e-5)) * ne * nHeII;

    // -- Free-free (Bremsstrahlung) cooling --
    // Lambda_ff = 1.42e-27 * g_ff(Z,T) * sqrt(T) * ne * n_ion * Z^2
    // van Hoof et al. (2014) Gaunt factor, computed per charge Z
    // gamma2 = Z^2 * Ry / (kT) where Ry = 2.1798723611035e-11 erg
    {
        constexpr double Ry_erg = 2.1798723611035e-11;
        double kT = kB_cgs * T;
        double sqrtT = sqrt(T);

        // Z=1: H+, He+
        double gamma2_Z1 = Ry_erg / kT;
        double gff_Z1 = totalGauntFactor(gamma2_Z1);
        cooling += 1.42e-27 * gff_Z1 * sqrtT * ne * (nHII + nHeII);

        // Z=2: He++
        double gamma2_Z2 = 4. * Ry_erg / kT;
        double gff_Z2 = totalGauntFactor(gamma2_Z2);
        cooling += 1.42e-27 * gff_Z2 * sqrtT * ne * 4. * nHeIII;

        // Metals: use per-charge Gaunt factor
        for (int e = 2; e < nElements; e++)
        {
            double nElem = nH * metalAbundances[e - 2];
            if (nElem < 1e-30) continue;
            int off = offset[e];
            int ns = effNS[e];
            for (int s = 1; s < ns; s++)  // skip neutral stage
            {
                double Z2 = double(s) * double(s);
                double gamma2 = Z2 * Ry_erg / kT;
                double gff = totalGauntFactor(gamma2);
                cooling += 1.42e-27 * gff * sqrtT * ne * nElem * ionFracs[off + s] * Z2;
            }
        }
    }

    // -- Metal line cooling (CHIANTI pre-tabulated) --
    // Lambda_metal = sum over elements and stages: n_elem * x_ion * Lambda_ion(T, ne)
    // where Lambda_ion is the per-ion cooling from the CHIANTI level population solver
    if (_hasCoolingTable)
    {
        for (int e = 2; e < nElements; e++)
        {
            double nElem = nH * metalAbundances[e - 2];
            if (nElem < 1e-30) continue;
            int off = offset[e];
            int ns = effNS[e];
            for (int s = 0; s < ns; s++)
            {
                int ionIdx = off + s;
                double xion = ionFracs[ionIdx];
                if (xion < 1e-30) continue;
                double LambdaIon = interpolateMetalCooling(ionIdx, T, ne);
                cooling += nElem * xion * LambdaIon;
            }
        }
    }

    // -- Charge exchange heating and cooling for metals --
    // H CX: recomb cooling (xi_HI) and ioniz heating (chi_HII) with dE = IP_s - IP_HI
    // He CX: recomb cooling (xi_HeI) and ioniz heating (chi_HeII) with dE = IP_s - IP_HeI
    // Kingdon & Ferland (1996) rates; net = cooling - heating
    {
        double IP_HI_erg = VernerCrossSections::HI_eV * eV_erg;
        double IP_HeI_erg = VernerCrossSections::HeI_eV * eV_erg;
        for (int e = 2; e < nElements; e++)
        {
            double nElem = nH * metalAbundances[e - 2];
            if (nElem < 1e-30) continue;
            int off = offset[e];
            int nxs = effNS[e] - 1;
            int voff = vernerOffset[e];

            for (int s = 0; s < nxs; s++)
            {
                int vi = voff + s;  // Verner index for stage s
                double IP_s_erg = VernerCrossSections::ionizationPotential(vi) * eV_erg;
                double x_next = ionFracs[off + s + 1];
                double x_this = ionFracs[off + s];

                // H charge exchange
                double dE_H = IP_s_erg - IP_HI_erg;
                if (x_next > 1e-30) cooling += dE_H * PhotoIonizationRates::xiHI(vi + 1, T) * nHI * x_next * nElem;
                if (x_this > 1e-30) cooling -= dE_H * PhotoIonizationRates::chiHII(vi, T) * nHII * x_this * nElem;

                // He charge exchange
                double dE_He = IP_s_erg - IP_HeI_erg;
                if (x_next > 1e-30) cooling += dE_He * PhotoIonizationRates::xiHeI(vi + 1, T) * nHeI * x_next * nElem;
                if (x_this > 1e-30) cooling -= dE_He * PhotoIonizationRates::chiHeII(vi, T) * nHeII * x_this * nElem;
            }
        }
    }

    return cooling;  // [erg s^-1 cm^-3]
}

//////////////////////////////////////////////////////////////////////

double PhotoIonizationSolver::solveIonizationBalance(double T, double ne, const double gamma[], double nH, double yHe,
                                                     const double metalAbundances[8], double ionFracs[],
                                                     const int effNS[]) const
{
    // iterate until electron density converges; start fully ionized (xe=1)
    double xe = 1.;
    ne = nH;
    double xePrev = xe;

    for (int iter = 0; iter < maxNeIterations; iter++)
    {
        // -- Hydrogen --
        double aHII = PhotoIonizationRates::alphaBHII(T);
        double bHI = PhotoIonizationRates::betaHI(T);
        double gHI = gamma[0];
        double cHI = (bHI + gHI / ne) / aHII;
        double xHI = 1. / (1. + cHI);
        double xHII = 1. - xHI;
        ionFracs[0] = xHI;
        ionFracs[1] = xHII;

        double HIe = xHI / xe;    // n_HI / n_e (per H density)
        double HIIe = xHII / xe;  // n_HII / n_e (per H density)

        // -- Helium --
        double aHeII = PhotoIonizationRates::alphaBHeII(T) + PhotoIonizationRates::alphaDRHeII(T);
        double bHeI = PhotoIonizationRates::betaHeI(T);
        double gHeI = gamma[1];
        double eHeII = PhotoIonizationRates::xiHIHeII(T);
        double aHeIII = PhotoIonizationRates::alphaBHeIII(T);
        double bHeII = PhotoIonizationRates::betaHeII(T);
        double gHeII = gamma[2];
        double eHeIII = PhotoIonizationRates::xiHIHeIII(T);

        double cHeI = safediv(bHeI + gHeI / ne, aHeII + eHeII * HIe);
        double cHeII = safediv(bHeII + gHeII / ne, aHeIII + eHeIII * HIe);
        double xHeI = 1. / (1. + cHeI * (1. + cHeII));
        double xHeII = xHeI * cHeI;
        double xHeIII = xHeII * cHeII;
        ionFracs[2] = xHeI;
        ionFracs[3] = xHeII;
        ionFracs[4] = xHeIII;

        double HeIe = xHeI * yHe / xe;    // n_HeI / n_e (per H density)
        double HeIIe = xHeII * yHe / xe;  // n_HeII / n_e (per H density)

        // -- Metals (Katz+ recursion) --
        xe = xHII + yHe * (xHeII + 2. * xHeIII);

        for (int e = 2; e < nElements; e++)
        {
            int off = offset[e];
            int ns = effNS[e];
            int nxs = ns - 1;  // number of ions with XS (= ns - 1, highest stage is sink)
            int voff = vernerOffset[e];

            // compute c_m for each ionization stage m (from highest to lowest)
            // c_m = (beta + chi_HII*HIIe + chi_HeII*HeIIe + gamma/ne) / (alpha + xi_HI*HIe + xi_HeI*HeIe)
            // Verner index voff+s corresponds to stage s (0=neutral); alpha/xi refer to stage s+1
            double c[20];  // max stages per element
            for (int s = nxs - 1; s >= 0; s--)
            {
                int vi = voff + s;  // Verner index for this ion stage
                double photoRate = gamma[vi];
                double betaRate = PhotoIonizationRates::beta(vi, T);
                double alphaRate = PhotoIonizationRates::alphaTotal(vi, T);  // recomb into stage s

                // charge exchange with H
                double xiHIRate = PhotoIonizationRates::xiHI(vi + 1, T);  // ion(s+1) + HI -> ion(s) + HII
                double chiHIIRate = PhotoIonizationRates::chiHII(vi, T);  // ion(s) + HII -> ion(s+1) + HI

                // charge exchange with He
                double xiHeIRate = PhotoIonizationRates::xiHeI(vi + 1, T);  // ion(s+1) + HeI -> ion(s) + HeII
                double chiHeIIRate = PhotoIonizationRates::chiHeII(vi, T);  // ion(s) + HeII -> ion(s+1) + HeI

                c[s] = safediv(betaRate + chiHIIRate * HIIe + chiHeIIRate * HeIIe + photoRate / ne,
                               alphaRate + xiHIRate * HIe + xiHeIRate * HeIe);
            }

            // recursion: x_0 = 1 / (1 + c_0 * (1 + c_1 * (1 + ... c_{N-1})))
            double denom = 1.;
            for (int s = nxs - 1; s >= 0; s--) denom = 1. + c[s] * denom;
            ionFracs[off] = 1. / denom;  // neutral fraction

            // forward fill; detect overflow (0 * Inf = NaN when c overflows)
            bool overflow = false;
            for (int s = 0; s < nxs; s++)
            {
                ionFracs[off + s + 1] = ionFracs[off + s] * c[s];
                if (std::isnan(ionFracs[off + s + 1]))
                {
                    overflow = true;
                    break;
                }
            }
            if (overflow)
            {
                // all c >> 1: gas is in highest ionization stage
                for (int s = 0; s < nxs; s++) ionFracs[off + s] = 0.;
                ionFracs[off + nxs] = 1.;
            }

            // accumulate electron contribution
            for (int s = 1; s < ns; s++) xe += metalAbundances[e - 2] * double(s) * ionFracs[off + s];
        }

        // stability: average with previous, enforce minimum electron fraction
        xe = 0.5 * (xe + xePrev);
        if (xe < minElectronFraction) xe = minElectronFraction;
        ne = xe * nH;

        // check convergence (absolute tolerance on xe)
        if (std::abs(xe - xePrev) < neTolerance) break;
        xePrev = xe;
    }

    return ne;
}

//////////////////////////////////////////////////////////////////////

PhotoIonizationSolver::CellResult PhotoIonizationSolver::solve(const Array& Jv, double nH, double yHe,
                                                               const double metalAbundances[8], double Tprior,
                                                               const double* ionFracsPrior, double TdampFactor) const
{
    CellResult result;

    // initialize ion fractions
    if (ionFracsPrior)
    {
        for (int i = 0; i < totalStages; i++) result.ionFracs[i] = ionFracsPrior[i];
    }
    else
    {
        initializeIonFractions(result.ionFracs);
    }

    // compute photoionization rates from J_nu
    double gamma[VernerCrossSections::numIons];
    computePhotoionizationRates(Jv, gamma);

    // initial electron density: always start fully ionized (xe=1)
    double ne = nH;

    // temperature bisection: find T where heating(T) = cooling(T)
    double Tlo = Tmin;
    double Thi = Tmax;

    // temperature damping: clamp bisection range to [Tprior/factor, Tprior*factor]
    if (TdampFactor > 1. && Tprior > Tmin)
    {
        Tlo = std::max(Tlo, Tprior / TdampFactor);
        Thi = std::min(Thi, Tprior * TdampFactor);
    }

    double T = std::max(Tlo, std::min(Tprior, Thi));
    double Tprev = T;

    for (int iter = 0; iter < maxTIterations; iter++)
    {
        // solve ionization balance at current T
        ne = solveIonizationBalance(T, ne, gamma, nH, yHe, metalAbundances, result.ionFracs, _nStagesEff);

        // compute heating and cooling
        double heat = computeHeatingRate(Jv, result.ionFracs, nH, yHe, metalAbundances, ne, _nStagesEff);
        double cool = computeCoolingRate(T, result.ionFracs, nH, yHe, metalAbundances, ne, _nStagesEff);

        // adaptive damping: tighten allowed T range as cell approaches equilibrium
        if (iter == 0 && TdampFactor > 1. && Tprior > Tmin && (heat + cool) > 0.)
        {
            double fractionalImbalance = std::abs(heat - cool) / (heat + cool);
            double adaptiveFactor = fractionalImbalance * (TdampFactor - 1.) + 1.;
            Tlo = std::max(Tlo, Tprior / adaptiveFactor);
            Thi = std::min(Thi, Tprior * adaptiveFactor);
        }

        // check convergence: relative T change (skip first two iterations to let bracket establish)
        if (iter > 1 && std::abs(T - Tprev) < Ttolerance * T)
        {
            result.Teq = T;
            result.ne = ne;
            result.heatingRate = heat;
            result.coolingRate = cool;
            return result;
        }

        // bisect: if heating > cooling, temperature needs to be higher (more cooling)
        // if cooling > heating, temperature needs to be lower (less cooling)
        if (heat > cool)
            Tlo = T;
        else
            Thi = T;

        // first step: jump to bracket edge to establish root-straddling bracket quickly
        Tprev = T;
        if (iter == 0)
            T = (heat > cool) ? Thi : Tlo;
        else
            // geometric bisection for wide brackets (>5 dex), linear otherwise
            T = (Thi > 1e5 * Tlo) ? std::sqrt(Tlo * Thi) : 0.5 * (Tlo + Thi);
    }

    // if not converged, return best estimate
    result.Teq = T;
    result.ne = ne;
    result.heatingRate = computeHeatingRate(Jv, result.ionFracs, nH, yHe, metalAbundances, ne, _nStagesEff);
    result.coolingRate = computeCoolingRate(T, result.ionFracs, nH, yHe, metalAbundances, ne, _nStagesEff);
    return result;
}

//////////////////////////////////////////////////////////////////////

PhotoIonizationSolver::CellResult PhotoIonizationSolver::solveIonizationAtFixedT(const Array& Jv, double nH, double yHe,
                                                                                 const double metalAbundances[8],
                                                                                 double T, const double* ionFracsPrior,
                                                                                 bool useEffectiveStages) const
{
    const int* effNS = useEffectiveStages ? _nStagesEff : nStages;

    CellResult result;

    // initialize ion fractions
    if (ionFracsPrior)
    {
        for (int i = 0; i < totalStages; i++) result.ionFracs[i] = ionFracsPrior[i];
    }
    else
    {
        initializeIonFractions(result.ionFracs);
    }

    // compute photoionization rates from J_lambda (included even in CIE mode)
    double gamma[VernerCrossSections::numIons];
    computePhotoionizationRates(Jv, gamma);

    // initial electron density: always start fully ionized (xe=1)
    double ne = nH;

    // solve ionization balance at fixed T (no bisection)
    ne = solveIonizationBalance(T, ne, gamma, nH, yHe, metalAbundances, result.ionFracs, effNS);

    // compute heating and cooling for diagnostics
    result.Teq = T;
    result.ne = ne;
    result.heatingRate = computeHeatingRate(Jv, result.ionFracs, nH, yHe, metalAbundances, ne, effNS);
    result.coolingRate = computeCoolingRate(T, result.ionFracs, nH, yHe, metalAbundances, ne, effNS);
    return result;
}

//////////////////////////////////////////////////////////////////////

double PhotoIonizationSolver::opacityAbs(double lambda, const double ionFracs[totalStages], double nH, double yHe,
                                         const double metalAbundances[8]) const
{
    double E_eV = Constants::h() * Constants::c() / (lambda * Constants::Qelectron());

    double kappa = 0.;

    // H
    kappa += nH * ionFracs[0] * VernerCrossSections::sigmaHI(E_eV);

    // He
    double nHe = nH * yHe;
    kappa += nHe * ionFracs[2] * VernerCrossSections::sigmaHeI(E_eV);
    kappa += nHe * ionFracs[3] * VernerCrossSections::sigmaHeII(E_eV);

    // metals
    for (int e = 2; e < nElements; e++)
    {
        double nElem = nH * metalAbundances[e - 2];
        int off = offset[e];
        int voff = vernerOffset[e];
        int nxs = nIonsXS[e];
        for (int s = 0; s < nxs; s++) kappa += nElem * ionFracs[off + s] * VernerCrossSections::sigma(voff + s, E_eV);
    }

    return kappa;  // [cm^-1]
}

//////////////////////////////////////////////////////////////////////
