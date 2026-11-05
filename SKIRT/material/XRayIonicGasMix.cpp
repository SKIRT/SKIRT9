/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "XRayIonicGasMix.hpp"
#include "AtomUtils.hpp"
#include "ComptonPhaseFunction.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "DipolePhaseFunction.hpp"
#include "ElectronScatteringHelper.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "LyUtils.hpp"
#include "MaterialState.hpp"
#include "NR.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"
#include "Range.hpp"
#include "StoredTable.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"
#include "VoigtProfile.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <limits>
#include <map>
#include <set>
#include <tuple>
#include <unordered_set>

////////////////////////////////////////////////////////////////////

namespace
{
    // ---- hardcoded configuration constants ----

    constexpr int numAtoms = 30;  // maximum atomic number used in this class

    // ---- common helper resources----

    // convert photon energy in eV to and from wavelength in m (same conversion in both directions)
    constexpr double wavelengthToFromEnergy(double x)
    {
        constexpr double front = Constants::h() * Constants::c() / Constants::Qelectron();
        return front / x;
    }

    // return true if the value is finite and represents an integer
    bool isFiniteInteger(double value)
    {
        return std::isfinite(value) && value == std::round(value)
               && value >= static_cast<double>(std::numeric_limits<int>::min())
               && value <= static_cast<double>(std::numeric_limits<int>::max());
    }

    // calculate the Voigt a parameter
    constexpr double voigtA(double lamGamma, double vth)
    {
        return lamGamma / (4. * M_PI * vth * M_SQRT2);
    }

    // Sample a spontaneous line in energy space, conditioned only on positive, finite energy.
    // The Lorentzian half-width is hbar*(Gamma_upper+Gamma_lower)/2 in eV. Do not truncate
    // the positive wings at the material cross-section limits or at an instrument's wavelength range.
    double sampleSpontaneousLineWavelength(double center, double energyHwhm, Random* random)
    {
        if (energyHwhm == 0.) return center;

        double energyCenter = wavelengthToFromEnergy(center);
        while (true)
        {
            double energy = energyCenter + energyHwhm * random->lorentz();
            if (energy > 0. && std::isfinite(energy)) return wavelengthToFromEnergy(energy);
        }
    }

    // wavelength range over which our cross sections may be nonzero
    constexpr Range nonZeroRange(wavelengthToFromEnergy(500e3), wavelengthToFromEnergy(4.3));

    // number of wavelengths per dex in high-resolution grid
    constexpr size_t numWavelengthsPerDex = 2500;

    // Return the dipole fraction of the angular redistribution matrix for an
    // electric-dipole transition J_l -> J_u.
    //
    // Following Hamilton (1947), the E1 redistribution matrix is written as a
    // weighted sum of the monopole and dipole terms. Here E1 and E2
    // denote the corresponding Hamilton coefficients for the given Delta J.
    // The returned value, E2 / (E1 + E2), is therefore the fraction of the
    // dipole component.
    constexpr double dipoleFractionFromJ(double lowerJ, double upperJ)
    {
        const double j = lowerJ;
        const double deltaJ = upperJ - lowerJ;

        double E1 = 0.0;
        double E2 = 0.0;

        if (deltaJ == 1)  // 1
        {
            E1 = 0.1 * 3.0 * j * (6.0 * j + 7.0) / ((j + 1.0) * (2.0 * j + 1.0));
            E2 = 0.1 * (2.0 * j + 5.0) * (j + 2.0) / ((j + 1.0) * (2.0 * j + 1.0));
        }
        else if (deltaJ == 0)  // 0
        {
            E1 = 0.1 * 3.0 * (2.0 * j * j + 2.0 * j + 1.0) / (j * (j + 1.0));
            E2 = 0.1 * (2.0 * j - 1.0) * (2.0 * j + 3.0) / (j * (j + 1.0));
        }
        else if (deltaJ == -1)  // -1
        {
            E1 = 0.1 * 3.0 * (j + 1.0) * (6.0 * j - 1.0) / (j * (2.0 * j + 1.0));
            E2 = 0.1 * (2.0 * j - 3.0) * (j - 1.0) / (j * (2.0 * j + 1.0));
        }
        else
        {
            E1 = 1.0;
            E2 = 0.0;
        }

        double sum = E1 + E2;
        if (sum <= 0.0) throw FATALERROR("Invalid monopole/dipole weights");

        double f = E2 / sum;
        if (f < 0.0 || f > 1.0) throw FATALERROR("Invalid dipole fraction");

        return f;
    }

    // resource data for photo-absorption
    struct PhotoAbsorbResource
    {
        PhotoAbsorbResource(const Array& a)
            : Z(a[0]), N(a[1]), n(a[2]), l(a[3]), Eth(a[4]), E0(a[5]), sigma0(a[6]), ya(a[7]), P(a[8]), yw(a[9])
        {}

        int ionIndex{-1};                    // index of the ion in the user-specified ion list
        int Z;                               // atomic number
        int N;                               // number of electrons
        int n;                               // principal quantum number of the shell
        int l;                               // orbital quantum number of the subshell
        double Eth;                          // subshell ionization threshold energy (eV)
        double E0;                           // fit parameter (eV)
        double sigma0;                       // fit parameter (Mb = 10^-22 m^2)
        double ya;                           // fit parameter (1)
        double P;                            // fit parameter (1)
        double yw;                           // fit parameter (1)
        static constexpr double Emax = 5e5;  // maximum energy for validity of the formula (eV)

        double Es{0.};        // thermal energy dispersion at the threshold energy (eV); set after construction
        double sigmamax{0.};  // cross section at the threshold energy plus 2*Es (m2); set after construction

        // return photo-absorption cross section in m2 for given energy in eV and cross section parameters,
        // without taking into account thermal dispersion
        double photoAbsorbSection(double E) const
        {
            if (E < Eth || E >= Emax) return 0.;

            double x = E / E0;
            double y = x;
            double xm1 = x - 1.;
            double Q = 5.5 + l - 0.5 * P;
            double F = (xm1 * xm1 + yw * yw) * std::pow(y, -Q) * std::pow(1. + std::sqrt(y / ya), -P);
            return 1e-22 * sigma0 * F;  // from Mb to m2
        }

        // return photo-absorption cross section in m2 for given energy in eV and cross section parameters,
        // approximating thermal dispersion by replacing the steep threshold transition by a sigmoid
        // error function with given parameters (dispersion and maximum value)
        double photoAbsorbThermalSection(double E) const
        {
            if (E <= Eth - 2. * Es) return 0.;
            if (E >= Eth + 2. * Es) return photoAbsorbSection(E);
            return sigmamax * (0.5 + 0.5 * std::erf((E - Eth) / Es));
        }
    };

    // resource data for fluorescence
    struct FluorescenceResource
    {
        FluorescenceResource(const Array& a) : Z(a[0]), N(a[1]), n(a[2]), l(a[3]), omega(a[4]), E(a[5]), W(a[6]) {}

        int paIndex{-1};  // index of the photo-absorption transition
        int Z;            // atomic number
        int N;            // number of electrons
        int n;            // principal quantum number of the shell with the hole
        int l;            // orbital quantum number of the subshell with the hole
        double omega;     // fluorescence yield (1)
        double E;         // (central) energy of the emitted photon (eV)
        double W;         // FWHM of the Lorentz shape for the emitted photon (eV), or zero
    };

    // resource data for line transitions
    struct LineResource
    {
        LineResource(const Array& a)
            : Z(a[0]), N(a[1]), lineIndex(a[2]), lowerLevelIndex(a.size() == 10 ? a[3] : 0),
              upperLevelIndex(a.size() == 10 ? a[4] : -1), lowerJ(a[a.size() == 10 ? 5 : 3]),
              upperJ(a[a.size() == 10 ? 6 : 4]), lam(a[a.size() == 10 ? 7 : 5]),
              A(a[a.size() == 10 ? 8 : 6]), lamGamma(a[a.size() == 10 ? 9 : 7])
        {}

        int ionIndex{-1};        // index of the corresponding ion in the user-specified ion list
        double scatterProb{-1};  // total probability that resonant absorption leads to re-emission
        double levelFraction{1.};  // population fraction of the entire lower level

        int Z;                // atomic number
        int N;                // number of electrons
        int lineIndex;        // persistent transition index
        int lowerLevelIndex;  // FAC lower index, or ground index 0 in the upstream format
        int upperLevelIndex;  // FAC upper index, or -1 when absent from the upstream format
        double lowerJ;        // lower total angular momentum
        double upperJ;        // upper total angular momentum
        double lam;           // wavelength (m)
        double A;             // Einstein A coefficient (1/s)
        double lamGamma;      // wavelength times natural decay rate (m/s)

        double section(double lambda, double vth) const
        {
            double a = voigtA(lamGamma, vth);
            double g = (2.0 * upperJ + 1.0) / (2.0 * lowerJ + 1.0);
            return LyUtils::section(lambda, lam, vth * M_SQRT2, A, a, g);
        }
    };

    // resource data for resonant-scattering branching probabilities
    struct BranchResource
    {
        BranchResource(const Array& a) : Z(a[0]), N(a[1]), upperIndex(a[2]), lowerIndex(a[3]), prob(a[4]) {}

        int ionIndex{-1};  // index of the corresponding ion in the user-specified ion list

        int Z;           // atomic number
        int N;           // number of electrons
        int upperIndex;  // index of the absorbed transition
        int lowerIndex;  // index of the emitted transition after branching
        double prob;     // branching probability from the upper transition to the lower transition
    };

    // resource data for FAC level-resolved photoionization events
    struct PhotoIonizationResource
    {
        PhotoIonizationResource(const Array& a)
        {
            for (int i = 0; i != 9; ++i)
                if (!isFiniteInteger(a[i]))
                    throw FATALERROR("FAC photoionization identifiers and quantum numbers must be finite integers");

            Z = static_cast<int>(a[0]);
            parentN = static_cast<int>(a[1]);
            targetN = static_cast<int>(a[2]);
            eventIndex = static_cast<int>(a[3]);
            parentBoundLevelIndex = static_cast<int>(a[4]);
            parentFreeLevelIndex = static_cast<int>(a[5]);
            targetLevelIndex = static_cast<int>(a[6]);
            n = static_cast<int>(a[7]);
            l = static_cast<int>(a[8]);
            Eth = a[9];
            E0 = a[10];
            sigma0 = a[11];
            ya = a[12];
            P = a[13];
            yw = a[14];
            y0 = a[15];
            y1 = a[16];

            if (Z <= 0 || Z > numAtoms || parentN <= 0 || parentN > Z || targetN != parentN - 1)
                throw FATALERROR("Invalid parent or target ion in FAC photoionization resource");
            if (eventIndex < 0 || parentBoundLevelIndex < 0 || parentFreeLevelIndex < 0 || targetLevelIndex < 0)
                throw FATALERROR("FAC photoionization event and level indices must be non-negative");
            if (n <= 0 || l < 0 || l >= n)
                throw FATALERROR("FAC photoionization subshell must satisfy n > 0 and 0 <= l < n");
            if (!std::isfinite(Eth) || !std::isfinite(E0) || !std::isfinite(sigma0) || !std::isfinite(ya)
                || !std::isfinite(P) || !std::isfinite(yw) || !std::isfinite(y0) || !std::isfinite(y1))
                throw FATALERROR("FAC photoionization fit parameters must be finite");
            if (Eth <= 0. || E0 <= 0. || sigma0 <= 0. || ya <= 0. || P <= 0. || yw < 0. || y0 < 0.
                || y1 <= 0.)
                throw FATALERROR("Invalid FAC photoionization fit parameters");
        }

        int ionIndex{-1};          // index of the parent ion in the user-specified ion list
        int iiParamIndex{-1};      // corresponding persistent II parameter, or -1 if no line branch exists
        double levelFraction{1.};  // population fraction of the parent bound level

        int Z;                       // atomic number
        int parentN;                 // electron number before photoionization
        int targetN;                 // electron number after photoionization
        int eventIndex;              // event index within the parent ion
        int parentBoundLevelIndex;   // FAC bound level of the parent ion
        int parentFreeLevelIndex;    // FAC continuum level of the parent ion
        int targetLevelIndex;        // matched FAC level of the target ion
        int n;                       // principal quantum number of the ionized subshell
        int l;                       // orbital quantum number of the ionized subshell
        double Eth;                  // ionization threshold energy (eV)
        double E0;                   // Verner-Yakovlev energy scale (eV)
        double sigma0;               // Verner-Yakovlev cross-section scale (Mb)
        double ya;                   // Verner-Yakovlev fit parameter (1)
        double P;                    // Verner-Yakovlev fit parameter (1)
        double yw;                   // Verner-Yakovlev fit parameter (1)
        double y0;                   // Verner et al. (1996) fit parameter (1)
        double y1;                   // Verner et al. (1996) fit parameter (1)

        static constexpr double Emax = 5e5;  // maximum energy for validity of the formula (eV)
        double Es;                          // thermal energy dispersion at the threshold (eV)
        double sigmamax;                    // cross section used to scale threshold smoothing (m2)

        // return the FAC event cross section in m2 without threshold thermal broadening
        double photoIonizationSection(double E) const
        {
            if (E < Eth || E >= Emax) return 0.;

            double x = E / E0 - y0;
            double y = std::hypot(x, y1);
            double xm1 = x - 1.;
            double logF = std::log(xm1 * xm1 + yw * yw) + (0.5 * P - 5.5 - l) * std::log(y)
                          - P * std::log1p(std::sqrt(y / ya));
            return std::exp(std::log(1e-22) + std::log(sigma0) + logF);  // from Mb to m2
        }

        // return the FAC event cross section in m2 with approximate threshold thermal broadening
        double photoIonizationThermalSection(double E) const
        {
            if (E <= Eth - 2. * Es) return 0.;
            if (E >= Eth + 2. * Es) return photoIonizationSection(E);
            return sigmamax * (0.5 + 0.5 * std::erf((E - Eth) / Es));
        }
    };

    // resource data for line branching following a FAC photoionization event
    struct InnerShellBranchResource
    {
        InnerShellBranchResource(const Array& a)
        {
            for (int i = 0; i != 5; ++i)
                if (!isFiniteInteger(a[i]))
                    throw FATALERROR("FAC II branching identifiers must be finite integers");

            Z = static_cast<int>(a[0]);
            parentN = static_cast<int>(a[1]);
            targetN = static_cast<int>(a[2]);
            eventIndex = static_cast<int>(a[3]);
            outputLineIndex = static_cast<int>(a[4]);
            prob = a[5];

            if (Z <= 0 || Z > numAtoms || parentN <= 0 || parentN > Z || targetN != parentN - 1)
                throw FATALERROR("Invalid parent or target ion in FAC II branching resource");
            if (eventIndex < 0 || outputLineIndex < 0)
                throw FATALERROR("FAC II event and output line indices must be non-negative");
            if (!std::isfinite(prob) || prob <= 0. || prob > 1.)
                throw FATALERROR("FAC II branching probability must be finite, larger than zero, and at most one");
        }

        int Z;                 // atomic number
        int parentN;           // electron number before photoionization
        int targetN;           // electron number after photoionization
        int eventIndex;        // referenced photoionization event
        int outputLineIndex;   // emitted transition in the target ion
        double prob;           // branching probability
    };

    // This function loads data from a resource file, only keeping the required data for the present ions.
    // The data is loaded into a vector of structs of type S that can be constructed from an array with C elements,
    // i.e. C columns.
    // The first two columns of the data must always be the ion parameters Z and N.
    // The function uses an set of unique ion indices to quickly check if an ion is present.
    // These indices must be calculated from the AtomUtils::ionIndex() function.
    template<class S, int C>
    vector<S> loadPresent(SimulationItem* item, string filename, string description,
                          const std::unordered_set<int>& ionSet)
    {
        vector<S> result;
        TextInFile infile(item, filename, description, true);
        infile.addColumn("Z");
        infile.addColumn("N");
        for (int i = 2; i < C; i++) infile.addColumn(string());

        Array row;
        while (infile.readRow(row))
        {
            if (!isFiniteInteger(row[0]) || !isFiniteInteger(row[1]))
                throw FATALERROR("Resource ion identifiers Z and N must be finite integers");

            int Z = row[0];
            int N = row[1];
            int index = AtomUtils::ionIndex(Z, N);

            // only keep ions that are present
            if (ionSet.count(index)) result.emplace_back(row);
        }
        return result;
    }

    // Accept the upstream eight-column format and its ten-column FAC extension.
    // A and the natural width are independent in both formats. Legacy lambda*A
    // resources must be converted explicitly instead of silently using A as Gamma.
    vector<LineResource> loadLines(SimulationItem* item, string filename, string description,
                                   const std::unordered_set<int>& ionSet)
    {
        std::ifstream stream(FilePaths::resource(filename));
        string text;
        while (std::getline(stream, text))
        {
            auto start = text.find_first_not_of(" \t\r");
            if (start == string::npos || text[start] == '#') continue;
            std::istringstream fields(text);
            double value;
            int count = 0;
            while (fields >> value) ++count;
            if (count == 8) return loadPresent<LineResource, 8>(item, filename, description, ionSet);
            if (count == 10) return loadPresent<LineResource, 10>(item, filename, description, ionSet);
            throw FATALERROR("Ionic_LN requires 8 upstream columns or 10 FAC columns with separate A and Gamma; "
                             "convert legacy 9-column lambda*A resources first");
        }
        throw FATALERROR("No numeric line records found in " + filename);
    }

    // return the resource filename for radiative recombination branching probabilities with a given number of electrons
    string branchRrFilename(int N)
    {
        return "Ionic_RR_N" + std::to_string(N) + ".stab";
    }
}

////////////////////////////////////////////////////////////////////

XRayIonicGasMix::XRayIonicGasMix(SimulationItem* parent, string ions, vector<double> abundances, double temperature,
                                 ElectronScattering electronScattering, bool resonantScattering, bool setup)
{
    _ions = ions;
    _abundances = abundances;
    _temperature = temperature;
    _electronScattering = electronScattering;
    _resonantScattering = resonantScattering;
    if (setup)
    {
        parent->addChild(this);
        this->setup();
    }
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMix::setupSelfBefore()
{
    MaterialMix::setupSelfBefore();

    auto config = find<Configuration>();

    // ------------ parse user properties ------------

    // parse all required ions from the ions property
    auto ionStrings = StringUtils::split(StringUtils::squeeze(ions()), ",");

    if (ionStrings.empty()) throw FATALERROR("No ions specified");
    if (ionStrings.size() != _abundances.size()) throw FATALERROR("Abundances must match amount of specified ions");

    // store all non-zero abundances and ions
    vector<double> abundances;
    abundances.reserve(_abundances.size());
    int numStrings = ionStrings.size();
    for (int i = 0; i < numStrings; i++)
    {
        if (_abundances[i] > 0.)
        {
            int Z, N;
            std::tie(Z, N) = AtomUtils::parseIon(ionStrings[i]);
            _ionParamv.push_back({Z, N});
            abundances.push_back(_abundances[i]);
        }
    }
    _numIons = _ionParamv.size();

    // store all unique N and ions
    std::set<int> usedN;
    std::unordered_set<int> ionSet;
    for (const auto& ion : _ionParamv)
    {
        usedN.insert(ion.N);
        ionSet.insert(AtomUtils::ionIndex(ion.Z, ion.N));
    }

    // create scattering helpers depending on the user-configured implementation type
    using namespace ElectronScatteringHelper;
    switch (electronScattering())
    {
        case ElectronScattering::None: _com = new NoScatteringHelper(this); break;
        case ElectronScattering::Free: _com = new FreeComptonHelper(this); break;
        case ElectronScattering::FreeWithPolarization: _com = new FreeComptonWithPolarizationHelper(this); break;
        case ElectronScattering::FreeBound: _com = new FreeBoundComptonHelper(this); break;
    }
    _numElec = electronScattering() == ElectronScattering::None ? 0 : _numIons;

    if (resonantScattering())
    {
        _dpf = new DipolePhaseFunction();
        _dpf->initialize(random(), true);
    }

    // ------------ load required resources (present ions only) ------------

    // These resources contain all the data needed during the setup, but are not needed afterwards.

    // photo-absorption data
    auto paResources = loadPresent<PhotoAbsorbResource, 10>(this, "Ionic_PA.txt", "photo-absorption data", ionSet);

    // fluorescence data
    auto fluoResources = loadPresent<FluorescenceResource, 7>(this, "Ionic_FL.txt", "fluorescence data", ionSet);
    vector<FluorescenceResource> heLikeIiResources;
    if (FilePaths::hasResource("Ionic_FL_II_N2.txt"))
        heLikeIiResources = loadPresent<FluorescenceResource, 7>(
            this, "Ionic_FL_II_N2.txt", "Li-like inner-shell-ionization to He-like z fixed-yield data", ionSet);

    // generic line data (recombination and resonant scattering)
    auto lineResources = loadLines(this, "Ionic_LN.txt", "resonant line data", ionSet);
    if (resonantScattering() && lineResources.empty())
        throw FATALERROR("Resonant scattering requires line data for at least one configured ion");

    // FAC photoionization events and their II branching data.
    vector<PhotoIonizationResource> piResources;
    if (FilePaths::hasResource("Ionic_PI.txt"))
        piResources = loadPresent<PhotoIonizationResource, 17>(
            this, "Ionic_PI.txt", "FAC photoionization event data", ionSet);
    vector<InnerShellBranchResource> iiBranchResources;
    if (FilePaths::hasResource("Ionic_BR_II.txt"))
        iiBranchResources = loadPresent<InnerShellBranchResource, 6>(
            this, "Ionic_BR_II.txt", "FAC II branching data", ionSet);

    // Load the line catalog for target ions populated by the selected parent-ion PI events.
    // This is separate from lineResources so that adding validation does not change RS, RR,
    // fluorescence, extinction, or scattering behavior.
    std::unordered_set<int> iiTargetIonSet;
    for (const auto& piRes : piResources) iiTargetIonSet.insert(AtomUtils::ionIndex(piRes.Z, piRes.targetN));
    auto iiLineResources =
        loadLines(this, "Ionic_LN.txt", "FAC II output line data", iiTargetIonSet);

    std::set<std::tuple<int, int, int, int>> piEventKeys;
    for (const auto& piRes : piResources)
    {
        auto key = std::make_tuple(piRes.Z, piRes.parentN, piRes.targetN, piRes.eventIndex);
        if (!piEventKeys.insert(key).second) throw FATALERROR("Duplicate FAC photoionization event key");
    }

    std::set<std::tuple<int, int, int>> iiLineKeys;
    for (const auto& lineRes : iiLineResources)
    {
        auto key = std::make_tuple(lineRes.Z, lineRes.N, lineRes.lineIndex);
        if (!iiLineKeys.insert(key).second) throw FATALERROR("Duplicate FAC II output line key");
    }

    std::set<std::tuple<int, int, int, int>> iiBranchKeys;
    std::map<std::tuple<int, int, int>, double> iiBranchSums;
    for (const auto& branchRes : iiBranchResources)
    {
        auto eventKey =
            std::make_tuple(branchRes.Z, branchRes.parentN, branchRes.targetN, branchRes.eventIndex);
        if (!piEventKeys.count(eventKey))
            throw FATALERROR("FAC II branching event is absent from the photoionization resource");

        auto lineKey = std::make_tuple(branchRes.Z, branchRes.targetN, branchRes.outputLineIndex);
        if (!iiLineKeys.count(lineKey)) throw FATALERROR("FAC II output line is absent from the line resource");

        auto branchKey =
            std::make_tuple(branchRes.Z, branchRes.parentN, branchRes.eventIndex, branchRes.outputLineIndex);
        if (!iiBranchKeys.insert(branchKey).second) throw FATALERROR("Duplicate FAC II branching key");

        auto sumKey = std::make_tuple(branchRes.Z, branchRes.parentN, branchRes.eventIndex);
        iiBranchSums[sumKey] += branchRes.prob;
    }
    for (const auto& item : iiBranchSums)
        if (item.second > 1. + 1e-6) throw FATALERROR("FAC II branching probabilities for an event sum to more than one");

    // Ions covered by FAC PI use the complete FAC event set instead of the legacy
    // shell-averaged PA/FL resources. This avoids counting the same photoionization twice.
    std::set<std::pair<int, int>> facPiIons;
    for (const auto& piRes : piResources) facPiIons.emplace(piRes.Z, piRes.parentN);
    paResources.erase(
        std::remove_if(paResources.begin(), paResources.end(),
                       [&facPiIons](const PhotoAbsorbResource& resource)
                       { return facPiIons.count(std::make_pair(resource.Z, resource.N)); }),
        paResources.end());
    fluoResources.erase(
        std::remove_if(fluoResources.begin(), fluoResources.end(),
                       [&facPiIons](const FluorescenceResource& resource)
                       { return facPiIons.count(std::make_pair(resource.Z, resource.N)); }),
        fluoResources.end());

    // Ground-configuration level populations used by resonant lines and FAC PI events.
    // If no file is provided, the ground level (FAC ilev=0) has unit population
    // and all other ground-configuration levels have zero population.
    if (resonantScattering() || !piResources.empty())
    {
        std::map<std::tuple<int, int, int>, double> levelFractions;

        if (!levelFractionsFilename().empty())
        {
            TextInFile infile(this, levelFractionsFilename(), "ground-configuration level population fractions");
            infile.addColumn("atomic number Z");
            infile.addColumn("electron number N");
            infile.addColumn("FAC level index");
            infile.addColumn("population fraction");

            std::map<std::pair<int, int>, double> fractionSums;
            Array row;
            while (infile.readRow(row))
            {
                if (!std::isfinite(row[0]) || !std::isfinite(row[1]) || !std::isfinite(row[2])
                    || row[0] != std::round(row[0]) || row[1] != std::round(row[1])
                    || row[2] != std::round(row[2]))
                    throw FATALERROR("Level population Z, N, and FAC level index must be finite integers");

                int Z = static_cast<int>(row[0]);
                int N = static_cast<int>(row[1]);
                int levelIndex = static_cast<int>(row[2]);
                double fraction = row[3];

                if (Z <= 0 || Z > numAtoms || N <= 0 || N > Z || levelIndex < 0)
                    throw FATALERROR("Invalid ion or FAC level index in level population file");
                if (!std::isfinite(fraction) || fraction < 0. || fraction > 1.)
                    throw FATALERROR("Level population fraction must be finite and between zero and one");

                auto key = std::make_tuple(Z, N, levelIndex);
                if (!levelFractions.emplace(key, fraction).second)
                    throw FATALERROR("Duplicate ion and FAC level index in level population file");
                fractionSums[std::make_pair(Z, N)] += fraction;
            }

            for (const auto& item : fractionSums)
            {
                if (std::abs(item.second - 1.) > 1e-8)
                    throw FATALERROR("Level population fractions for an ion must sum to one");
            }
        }
        else
        {
            if (resonantScattering())
            {
                for (const auto& lineRes : lineResources)
                {
                    auto key = std::make_tuple(lineRes.Z, lineRes.N, lineRes.lowerLevelIndex);
                    levelFractions[key] = lineRes.lowerLevelIndex == 0 ? 1. : 0.;
                }
            }
            for (const auto& piRes : piResources)
            {
                auto key = std::make_tuple(piRes.Z, piRes.parentN, piRes.parentBoundLevelIndex);
                levelFractions[key] = piRes.parentBoundLevelIndex == 0 ? 1. : 0.;
            }
        }

        if (resonantScattering())
        {
            for (auto& lineRes : lineResources)
            {
                auto fraction =
                    levelFractions.find(std::make_tuple(lineRes.Z, lineRes.N, lineRes.lowerLevelIndex));
                if (fraction == levelFractions.end())
                    throw FATALERROR("Missing lower-level population for a resonant line");
                lineRes.levelFraction = fraction->second;
            }
        }
        for (auto& piRes : piResources)
        {
            auto fraction =
                levelFractions.find(std::make_tuple(piRes.Z, piRes.parentN, piRes.parentBoundLevelIndex));
            if (fraction == levelFractions.end())
                throw FATALERROR("Missing parent-level population for a FAC photoionization event");
            piRes.levelFraction = fraction->second;
        }
    }

    // branching probability data
    vector<BranchResource> branchResources;
    if (resonantScattering())
        branchResources =
            loadPresent<BranchResource, 5>(this, "Ionic_BR.txt", "resonant branching probabilities", ionSet);

    // yields for recombination for each value of N
    std::map<int, StoredTable<3>> recoResources;
    for (int N : usedN)
    {
        string rrFilename = branchRrFilename(N);
        if (FilePaths::hasResource(rrFilename))
            recoResources.emplace(N, StoredTable<3>(this, rrFilename, "Z(1),Index(1),T(K)", "Y(1)"));
    }

    if (recoResources.size() == 0) find<Log>()->warning("No recombination resources found");

    // ------------ preprocess resources ------------

    // Recombination lines can be modelled as scattering following an inner-shell photo-absorption (n,l)=(1,0).
    // This ignores the cascade and only models the transition back to the inner shell.
    // We can thus model it the same as fluorescence and will add it to the fluorescence resources.
    for (const auto& lineRes : lineResources)
    {
        auto recoIt = recoResources.find(lineRes.N);
        if (recoIt == recoResources.end()) continue;

        double Z = lineRes.Z;
        double N = lineRes.N;
        double lineIndex = lineRes.lineIndex;
        double E = wavelengthToFromEnergy(lineRes.lam);
        double omega = recoIt->second(Z, lineIndex, temperature());

        if (omega == 0) continue;

        Array params = {Z, N, 1., 0., omega, E, 0.};
        fluoResources.emplace_back(params);
    }

    // Li-like inner-shell ionization -> He-like z.
    // Keep these SPEX-based fixed-yield lines separate from the FAC Li-like line catalog.
    fluoResources.insert(fluoResources.end(), heLikeIiResources.begin(), heLikeIiResources.end());

    // ------------ postprocess resources ------------

    int numPa = paResources.size();
    _numFluo = fluoResources.size();
    _numLine = lineResources.size();

    // store the cross reference indices for pa, fluo, and res
    for (int i = 0; i < _numIons; i++)
    {
        auto& ion = _ionParamv[i];

        // reference ion in each pa
        for (int p = 0; p < numPa; p++)
        {
            auto& paRes = paResources[p];
            if (paRes.Z == ion.Z && paRes.N == ion.N)
            {
                paRes.ionIndex = i;

                // reference pa in each fluo
                for (auto& fluoRes : fluoResources)
                {
                    if (fluoRes.Z == paRes.Z && fluoRes.N == paRes.N && fluoRes.n == paRes.n && fluoRes.l == paRes.l)
                        fluoRes.paIndex = p;
                }
            }
        }

        // reference parent ion in each FAC PI event
        for (auto& piRes : piResources)
        {
            if (piRes.Z == ion.Z && piRes.parentN == ion.N) piRes.ionIndex = i;
        }

        // reference ion in each res
        if (resonantScattering())
        {
            for (auto& lineRes : lineResources)
            {
                if (lineRes.Z == ion.Z && lineRes.N == ion.N) lineRes.ionIndex = i;
            }
            for (auto& braRes : branchResources)
            {
                if (braRes.Z == ion.Z && braRes.N == ion.N) braRes.ionIndex = i;
            }
        }
    }

    // calculate the persistent thermal velocities
    // doesn't actually use any resources, but stores this for all 30 atomic numbers
    _vthermv.resize(numAtoms, 0.);
    for (int Z = 1; Z <= numAtoms; Z++) _vthermv[Z - 1] = sqrt(Constants::k() * temperature() / AtomUtils::mass(Z));

    // Photo-absorption
    // calculate the parameters for the sigmoid function approximating the convolution with a Gaussian
    // at the threshold energy for each cross section record, and store the result into a temporary vector;
    // the information includes the thermal energy dispersion at the threshold energy and
    // the intrinsic cross section at the threshold energy plus twice this energy dispersion
    for (auto& paRes : paResources)
    {
        const auto& ion = _ionParamv[paRes.ionIndex];
        paRes.Es = paRes.Eth * vtherm(ion.Z) / Constants::c();
        paRes.sigmamax = paRes.photoAbsorbSection(paRes.Eth + 2. * paRes.Es);
    }
    for (auto& piRes : piResources)
    {
        if (piRes.ionIndex < 0) throw FATALERROR("FAC photoionization event has no configured parent ion");
        const auto& ion = _ionParamv[piRes.ionIndex];
        piRes.Es = piRes.Eth * vtherm(ion.Z) / Constants::c();
        piRes.sigmamax = piRes.photoIonizationSection(piRes.Eth + 2. * piRes.Es);
    }

    // Resonant scattering
    // Calculate the total branching probability for each resonant transition.
    // This is the probability that a new photon will be emitted after a 'resonant absorption' event.
    // This value can be lower than 1 because low-energy photons are ignored and thus not re-emitted.
    if (resonantScattering())
    {
        for (auto& lineRes : lineResources)
        {
            // total branching probability
            lineRes.scatterProb = 0.;
            for (const auto& braRes : branchResources)
            {
                if (braRes.ionIndex == lineRes.ionIndex && braRes.upperIndex == lineRes.lineIndex)
                    lineRes.scatterProb += braRes.prob;

                if (lineRes.scatterProb == 1.) break;  // speed up since branching matrix has a lot of 1s and 0s
            }
        }
    }

    // ------------ calculate/store persistent data ------------

    // The persistent data is the data that is needed beyond the setup (scattering)
    // No changes should be made to the resource vectors after this point.
    // Their ordering must remain consistent with the corresponding persistent parameter vectors.

    // Fluorescence
    // Store the vtherm, wavelength, and width of each fluorescence transition.
    // These are needed during scattering after a photo-absorption event.
    _fluorescenceParamv.resize(_numFluo);
    for (int f = 0; f != _numFluo; f++)
    {
        const auto& fluoRes = fluoResources[f];
        auto& fluo = _fluorescenceParamv[f];

        fluo.vth = vtherm(fluoRes.Z);
        fluo.lambda = wavelengthToFromEnergy(fluoRes.E);
        fluo.width = fluoRes.W / 2.;  // convert from FWHM to HWHM
    }

    // Common output-line catalog for RS branching and FAC II. Start with all configured-ion
    // lines used by RS, RR, and legacy processing, then append target-ion lines required only
    // by FAC II. Duplicate (Z,N,line_index) entries refer to the same catalog position.
    vector<LineResource> outputLineResources = lineResources;
    std::set<std::tuple<int, int, int>> outputLineKeys;
    for (const auto& lineRes : outputLineResources)
    {
        auto key = std::make_tuple(lineRes.Z, lineRes.N, lineRes.lineIndex);
        if (!outputLineKeys.insert(key).second) throw FATALERROR("Duplicate line in common output-line catalog");
    }
    for (const auto& lineRes : iiLineResources)
    {
        auto key = std::make_tuple(lineRes.Z, lineRes.N, lineRes.lineIndex);
        if (outputLineKeys.insert(key).second) outputLineResources.push_back(lineRes);
    }

    std::map<std::tuple<int, int, int>, int> linePosition;
    _lineParamv.resize(outputLineResources.size());
    for (size_t l = 0; l != outputLineResources.size(); ++l)
    {
        const auto& lineRes = outputLineResources[l];
        auto& line = _lineParamv[l];

        line.Z = lineRes.Z;
        line.N = lineRes.N;
        line.lineIndex = lineRes.lineIndex;
        line.vth = vtherm(lineRes.Z);
        line.lambda = lineRes.lam;
        line.a = voigtA(lineRes.lamGamma, line.vth);
        // E0*lambda*Gamma/(4*pi*c) = hbar*Gamma/2, independent of the thermal velocity.
        line.naturalEnergyHwhm =
            wavelengthToFromEnergy(lineRes.lam) * lineRes.lamGamma / (4. * M_PI * Constants::c());
        if (!std::isfinite(line.naturalEnergyHwhm) || line.naturalEnergyHwhm < 0.)
            throw FATALERROR("Invalid natural energy width in output-line catalog");

        auto key = std::make_tuple(line.Z, line.N, line.lineIndex);
        if (!linePosition.emplace(key, static_cast<int>(l)).second)
            throw FATALERROR("Duplicate position in common output-line catalog");
    }

    // FAC inner-shell ionization
    // Store each PI event with line re-emission and its explicit output-line positions
    // in the common catalog. The cumulative distribution is conditional on line
    // re-emission; scatterProb retains the absolute total branching probability.
    size_t storedIiBranches = 0;
    for (auto& piRes : piResources)
    {
        vector<double> branchingProbabilities;
        vector<int> branchLinePositions;
        std::set<int> outputPositions;

        for (const auto& branchRes : iiBranchResources)
        {
            if (branchRes.Z != piRes.Z || branchRes.parentN != piRes.parentN
                || branchRes.targetN != piRes.targetN || branchRes.eventIndex != piRes.eventIndex)
                continue;

            auto output =
                linePosition.find(std::make_tuple(branchRes.Z, branchRes.targetN, branchRes.outputLineIndex));
            if (output == linePosition.end())
                throw FATALERROR("FAC II output line is absent from the common output-line catalog");
            if (!outputPositions.insert(output->second).second)
                throw FATALERROR("Duplicate FAC II output line for a photoionization event");

            branchLinePositions.push_back(output->second);
            branchingProbabilities.push_back(branchRes.prob);
        }

        if (branchingProbabilities.empty()) continue;

        piRes.iiParamIndex = static_cast<int>(_innerShellIonizationParamv.size());
        _innerShellIonizationParamv.emplace_back();
        auto& ii = _innerShellIonizationParamv.back();
        ii.Z = piRes.Z;
        ii.parentN = piRes.parentN;
        ii.eventIndex = piRes.eventIndex;
        ii.branchLinePositionv = std::move(branchLinePositions);

        Array probabilities(branchingProbabilities.size());
        for (size_t b = 0; b != branchingProbabilities.size(); ++b)
            probabilities[b] = branchingProbabilities[b];
        ii.scatterProb = NR::cdf(ii.cumBranchingv, probabilities);
        if (ii.scatterProb <= 0. || ii.scatterProb > 1. + 1e-6)
            throw FATALERROR("Invalid total FAC II branching probability");

        storedIiBranches += branchingProbabilities.size();
    }
    if (_innerShellIonizationParamv.size() != iiBranchSums.size())
        throw FATALERROR("FAC II event count is inconsistent after persistent branching setup");
    if (storedIiBranches != iiBranchResources.size())
        throw FATALERROR("FAC II branch count is inconsistent after persistent branching setup");
    _numII = static_cast<int>(_innerShellIonizationParamv.size());

    // Resonant scattering
    // Store each resonant input line and its explicit output-line positions in the common catalog.
    if (resonantScattering())
    {
        _resonantParamv.resize(_numLine);

        for (int r = 0; r != _numLine; r++)
        {
            const auto& lineRes = lineResources[r];
            auto& res = _resonantParamv[r];

            res.ionIndex = lineRes.ionIndex;
            res.lineIndex = lineRes.lineIndex;
            res.dipoleFraction = dipoleFractionFromJ(lineRes.lowerJ, lineRes.upperJ);

            auto input = linePosition.find(std::make_tuple(lineRes.Z, lineRes.N, lineRes.lineIndex));
            if (input == linePosition.end())
                throw FATALERROR("Resonant input line is absent from the common output-line catalog");
            res.inputLinePosition = input->second;

            // Store explicit output-line positions and their conditional CDF.
            // This supports sparse indices as well as branching in either index direction.
            vector<double> branchingProbabilities;
            std::set<int> outputPositions;
            for (const auto& braRes : branchResources)
            {
                if (braRes.ionIndex != res.ionIndex || braRes.upperIndex != res.lineIndex) continue;

                if (braRes.prob <= 0. || braRes.prob > 1.)
                    throw FATALERROR("Resonant branching probability must be larger than zero and at most one");

                auto output = linePosition.find(std::make_tuple(braRes.Z, braRes.N, braRes.lowerIndex));
                if (output == linePosition.end())
                    throw FATALERROR("Resonant branching output line is absent from the common output-line catalog");
                if (!outputPositions.insert(output->second).second)
                    throw FATALERROR("Duplicate resonant branching output for an input line");

                res.branchLinePositionv.push_back(output->second);
                branchingProbabilities.push_back(braRes.prob);
            }

            if (!branchingProbabilities.empty())
            {
                Array scatProb(branchingProbabilities.size());
                for (size_t b = 0; b != branchingProbabilities.size(); ++b) scatProb[b] = branchingProbabilities[b];
                double totalProbability = NR::cdf(res.cumBranchingv, scatProb);
                if (totalProbability > 1. + 1e-6)
                    throw FATALERROR("Resonant branching probabilities for an input line sum to more than one");
            }
        }
    }

    // ------------ wavelength grid ------------

    // construct a wavelength grid for sampling cross sections containing a merged set of grid points
    // in the relevant wavelength range (intersection of simulation range and nonzero range):
    //  - a fine grid in log space that provides sufficient resolution for most applications
    //  - all specific wavelengths mentioned in the configuration of the simulation (grids, normalizations, ...)
    //    ensuring that the cross sections are calculated at exactly these wavelengths
    //  - 7 extra wavelength points around the threshold energies for all transitions,
    //    placed at -2, -4/3, -2/3, 0, 2/3, 4/3, 2 times the thermal energy dispersion

    // get the relevant range (intersection of simulation range and nonzero range)
    Range range = config->simulationWavelengthRange();
    range.intersect(nonZeroRange);

    // we first gather all the wavelength points, in arbitrary order, and then sort them
    vector<double> lambdav;
    lambdav.reserve(5 * numWavelengthsPerDex);

    // add a fine grid in log space;
    // use integer multiples as logarithmic grid points so that the grid is stable for changing wavelength ranges
    constexpr double numPerDex = numWavelengthsPerDex;  // converted to double to avoid casting
    int minLambdaSerial = std::floor(numPerDex * log10(range.min()));
    int maxLambdaSerial = std::ceil(numPerDex * log10(range.max()));
    for (int k = minLambdaSerial; k <= maxLambdaSerial; k++) lambdav.push_back(pow(10., k / numPerDex));

    // add the wavelengths mentioned in the configuration of the simulation
    for (double lambda : config->simulationWavelengths())
        if (range.contains(lambda)) lambdav.push_back(lambda);

    // add wavelength points around the threshold energies for all transitions
    for (const auto& paRes : paResources)
    {
        double Es = paRes.Es;
        for (double delta : {-2., -4. / 3., -2. / 3., 0., 2. / 3., 4. / 3., 2.})
        {
            double lambda = wavelengthToFromEnergy(paRes.Eth + delta * Es);
            if (range.contains(lambda)) lambdav.push_back(lambda);
        }
    }
    for (const auto& piRes : piResources)
    {
        double Es = piRes.Es;
        for (double delta : {-2., -4. / 3., -2. / 3., 0., 2. / 3., 4. / 3., 2.})
        {
            double lambda = wavelengthToFromEnergy(piRes.Eth + delta * Es);
            if (range.contains(lambda)) lambdav.push_back(lambda);
        }
    }

    // add the fluorescence emission wavelengths
    for (const auto& fluoRes : fluoResources)
    {
        double lambda = wavelengthToFromEnergy(fluoRes.E);
        if (range.contains(lambda)) lambdav.push_back(lambda);
    }

    // add the outer wavelengths of our nonzero range, plus an extra just outside of that range,
    // so that there are always at least three points and thus two bins in the grid
    lambdav.push_back(nonZeroRange.min());
    lambdav.push_back(nonZeroRange.max());
    lambdav.push_back(nonZeroRange.max() * 1.000001);  // this wavelength point is never actually used

    // sort the wavelengths and remove duplicates
    NR::unique(lambdav);
    int numLambda = lambdav.size();

    // derive a wavelength grid that will be used for converting a wavelength to an index in the above array;
    // the grid points are shifted to the left of the actual sample points to approximate rounding
    _lambdav.resize(numLambda);
    _lambdav[0] = lambdav[0];
    for (int ell = 1; ell != numLambda; ell++)
    {
        _lambdav[ell] = sqrt(lambdav[ell] * lambdav[ell - 1]);
    }

    // ------------ extinction ------------

    // calculate the extinction cross section at every wavelength; to guarantee that the cross section is zero
    // for wavelengths outside our range, leave the values for the three outer wavelength points at zero
    _sigmaextv.resize(numLambda);
    for (int ell = 1; ell < numLambda - 2; ell++)
    {
        double lambda = lambdav[ell];
        double sigma = 0.;

        // electron scattering
        for (int i = 0; i < _numElec; i++)
        {
            const auto& ion = _ionParamv[i];
            sigma += _com->sectionSca(lambda, ion.Z, ion.N) * abundances[i];
        }

        // photo-absorption and fluorescence
        for (const auto& paRes : paResources)
        {
            double E = wavelengthToFromEnergy(lambda);
            sigma += paRes.photoAbsorbThermalSection(E) * abundances[paRes.ionIndex];
        }
        for (const auto& piRes : piResources)
        {
            double E = wavelengthToFromEnergy(lambda);
            sigma += piRes.photoIonizationThermalSection(E) * abundances[piRes.ionIndex]
                     * piRes.levelFraction;
        }

        // resonant scattering
        if (resonantScattering())
        {
            for (const auto& lineRes : lineResources)
            {
                double vth = vtherm(lineRes.Z);
                sigma += lineRes.section(lambda, vth) * abundances[lineRes.ionIndex] * lineRes.levelFraction;
            }
        }

        _sigmaextv[ell] = sigma;
    }

    // ------------ scattering ------------

    // make room for the scattering cross section and the cumulative fluorescence/scattering probabilities
    _sigmascav.resize(numLambda);
    _cumsigmascavv.resize(numLambda, 0);

    // provide temporary array for the non-normalized fluorescence/scattering contributions (at the current wavelength)
    int numInteractions = _numElec + _numFluo + _numLine + _numII;
    Array sigmas(numInteractions);

    // calculate the above for every wavelength; as before, leave the values for the outer wavelength points at zero
    for (int ell = 1; ell < numLambda - 2; ell++)
    {
        double lambda = lambdav[ell];
        double E = wavelengthToFromEnergy(lambda);

        // electron scattering
        for (int i = 0; i < _numElec; i++)
        {
            const auto& ion = _ionParamv[i];

            sigmas[i] = _com->sectionSca(lambda, ion.Z, ion.N) * abundances[i];
        }

        // fluorescence: iterate over both cross section and fluorescence parameter sets in sync
        for (int f = 0; f < _numFluo; f++)
        {
            const auto& fluoRes = fluoResources[f];
            const auto& paRes = paResources[fluoRes.paIndex];

            double sigma = paRes.photoAbsorbThermalSection(E) * abundances[paRes.ionIndex] * fluoRes.omega;
            sigmas[_numElec + f] = sigma;
        }

        // resonant scattering
        if (resonantScattering())
        {
            for (int r = 0; r < _numLine; r++)
            {
                const auto& lineRes = lineResources[r];

                double vth = vtherm(lineRes.Z);
                double sigma = lineRes.section(lambda, vth) * abundances[lineRes.ionIndex]
                               * lineRes.levelFraction * lineRes.scatterProb;
                sigmas[_numElec + _numFluo + r] = sigma;
            }
        }

        // FAC inner-shell ionization
        for (const auto& piRes : piResources)
        {
            if (piRes.iiParamIndex < 0) continue;

            const auto& ii = _innerShellIonizationParamv[piRes.iiParamIndex];
            double sigma = piRes.photoIonizationThermalSection(E) * abundances[piRes.ionIndex]
                           * piRes.levelFraction * ii.scatterProb;
            sigmas[_numElec + _numFluo + _numLine + piRes.iiParamIndex] = sigma;
        }

        // determine the normalized cumulative probability distribution and the cross section
        _sigmascav[ell] = NR::cdf(_cumsigmascavv[ell], sigmas);
    }
}

////////////////////////////////////////////////////////////////////

XRayIonicGasMix::~XRayIonicGasMix()
{
    delete _com;
    delete _dpf;
}

////////////////////////////////////////////////////////////////////

int XRayIonicGasMix::indexForLambda(double lambda) const
{
    return NR::locateFail(_lambdav, lambda);
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::vtherm(int Z) const
{
    return _vthermv[Z - 1];
}

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType XRayIonicGasMix::materialType() const
{
    return MaterialMix::MaterialType::Gas;
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::hasPolarizedScattering() const
{
    return electronScattering() == ElectronScattering::FreeWithPolarization || resonantScattering();
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::hasResonantScattering() const
{
    return resonantScattering();
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::hasScatteringDispersion() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::scatteringEmulatesSecondaryEmission() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> XRayIonicGasMix::specificStateVariableInfo() const
{
    return vector<StateVariable>{StateVariable::numberDensity()};
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::mass() const
{
    return Constants::Mproton();
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::sectionAbs(double lambda) const
{
    int index = indexForLambda(lambda);
    if (index < 0) return 0.;
    return _sigmaextv[index] - _sigmascav[index];
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::sectionSca(double lambda) const
{
    int index = indexForLambda(lambda);
    if (index < 0) return 0.;
    return _sigmascav[index];
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::sectionExt(double lambda) const
{
    int index = indexForLambda(lambda);
    if (index < 0) return 0.;
    return _sigmaextv[index];
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double number = state->numberDensity();
    return number > 0. ? sectionAbs(lambda) * number : 0.;
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::opacitySca(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double number = state->numberDensity();
    return number > 0. ? sectionSca(lambda) * number : 0.;
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double number = state->numberDensity();
    return number > 0. ? sectionExt(lambda) * number : 0.;
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMix::setScatteringInfoIfNeeded(PhotonPacket* pp, const MaterialState* state, const double lambda) const
{
    auto scatinfo = pp->getScatteringInfo();

    if (!scatinfo->valid)
    {
        scatinfo->valid = true;

        // indexForLambda should never be able to fail here, check needed anyway?
        scatinfo->species = NR::locateClip(_cumsigmascavv[indexForLambda(lambda)], random()->uniform());

        // Compton scattering
        if (scatinfo->species < _numElec)
        {
            int i = scatinfo->species;
            const auto& ion = _ionParamv[i];
            scatinfo->velocity = vtherm(ion.Z) * random()->maxwell();
        }
        // Fluorescent emission (scattering)
        else if (scatinfo->species < _numElec + _numFluo)
        {
            int f = scatinfo->species - _numElec;
            const auto& fluo = _fluorescenceParamv[f];
            double lambda = fluo.lambda;
            double width = fluo.width;

            scatinfo->velocity = fluo.vth * random()->maxwell();

            if (width == 0.)
            {
                // for a zero-width line, simply copy the central wavelength
                scatinfo->lambda = lambda;
            }
            else
            {
                // otherwise sample a wavelength from the Lorentz line shape in energy space;
                // the tails of the Lorentz distribition are very long, occasionaly resulting in negative energies;
                // therefore we loop until the sampled wavelength is meaningful
                while (true)
                {
                    double center = wavelengthToFromEnergy(lambda);
                    scatinfo->lambda = wavelengthToFromEnergy(center + width * random()->lorentz());
                    if (nonZeroRange.contains(scatinfo->lambda)) break;
                }
            }
        }
        // Resonant scattering
        else if (scatinfo->species < _numElec + _numFluo + _numLine)
        {
            int r = scatinfo->species - _numElec - _numFluo;
            const auto& res = _resonantParamv[r];
            const auto& inputLine = _lineParamv[res.inputLinePosition];

            int branch = NR::locateFail(res.cumBranchingv, random()->uniform());

            if (branch == -1) throw FATALERROR("Sampling from resonant branching probability has failed");

            int lr = res.branchLinePositionv[branch];
            if (lr < 0 || lr >= static_cast<int>(_lineParamv.size()))
                throw FATALERROR("Resonant branching line position is out of range");

            // if coherent (no branching)
            if (lr == res.inputLinePosition)
            {
                scatinfo->lambda = 0.;  // inform scatter is coherent

                // determine whether this scattering event uses the dipole phase function
                scatinfo->dipole = random()->uniform() < res.dipoleFraction;

                // sample an atom velocity from the Voigt profile (the atom-frame frequency it also
                // returns is not needed here since the dipole choice above already used dipoleFraction)
                scatinfo->velocity =
                    LyUtils::sampleAtomVelocity(lambda, inputLine.lambda, M_SQRT2 * inputLine.vth, inputLine.a,
                                                temperature(), state->numberDensity(), pp->direction(), config(), random())
                        .first;
            }
            // if incoherent (branching)
            else
            {
                const auto& outputLine = _lineParamv[lr];

                // a cascade erases the incoming direction and the atom-frame frequency, but not the velocity
                // of the absorbing atom, which also emits the photon; so we sample this velocity for the
                // absorbing (upper) transition as for coherent scattering, i.e. with the component along the
                // incoming direction drawn from the Voigt profile and the perpendicular components drawn from
                // the Maxwell distribution, but without the acceleration scheme, which is intended for coherent
                // scattering in the line core; the photon is then emitted isotropically at the wavelength of
                // the lower transition, Doppler shifted by this velocity
                Direction kin = pp->direction();
                double vD = M_SQRT2 * inputLine.vth;
                double x = (inputLine.lambda - lambda) / lambda * Constants::c() / vD;
                double upar = VoigtProfile::sample(inputLine.a, x, random());
                Vec velocity = inputLine.vth * random()->maxwell();
                velocity += (upar * vD - Vec::dot(velocity, kin)) * kin;

                scatinfo->lambda = outputLine.lambda;
                scatinfo->dipole = false;
                scatinfo->velocity = velocity;
            }
        }
        // Line emission following FAC inner-shell ionization
        else
        {
            int iiIndex = scatinfo->species - _numElec - _numFluo - _numLine;
            if (iiIndex < 0 || iiIndex >= _numII) throw FATALERROR("FAC II scattering event index is out of range");
            const auto& ii = _innerShellIonizationParamv[iiIndex];

            int branch = NR::locateFail(ii.cumBranchingv, random()->uniform());
            if (branch == -1) throw FATALERROR("Sampling from FAC II branching probability has failed");

            int linePosition = ii.branchLinePositionv[branch];
            if (linePosition < 0 || linePosition >= static_cast<int>(_lineParamv.size()))
                throw FATALERROR("FAC II branching line position is out of range");
            const auto& outputLine = _lineParamv[linePosition];

            // Direct and cascade-fed II lines use the width of the selected emitting transition.
            // Cache the atom-frame wavelength here so all peel-offs and the actual packet share it.
            scatinfo->lambda =
                sampleSpontaneousLineWavelength(outputLine.lambda, outputLine.naturalEnergyHwhm, random());
            scatinfo->velocity = outputLine.vth * random()->maxwell();
        }
    }
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs,
                                        Direction bfky, const MaterialState* state, const PhotonPacket* pp) const
{
    // draw a random scattering channel and atom velocity, unless a previous peel-off stored this already
    setScatteringInfoIfNeeded(const_cast<PhotonPacket*>(pp), state, lambda);
    auto scatinfo = const_cast<PhotonPacket*>(pp)->getScatteringInfo();

    // Compton scattering in electron rest frame; with support for polarization if enabled
    if (scatinfo->species < _numElec)
    {
        int i = scatinfo->species;
        const auto& ion = _ionParamv[i];
        // transform the wavelength into the rest frame of the electron
        lambda = PhotonPacket::shiftedReceptionWavelength(lambda, pp->direction(), scatinfo->velocity);
        _com->peeloffScattering(I, Q, U, V, lambda, ion.Z, ion.N, pp->direction(), bfkobs, bfky, pp);
        lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfkobs, scatinfo->velocity);
        return false;
    }

    // fluorescence
    else if (scatinfo->species < _numElec + _numFluo)
    {
        // unpolarized isotropic emission; the bias weight is trivially 1 and there is no contribution to Q, U, V
        I = 1.;

        // update the photon packet wavelength to the (possibly sampled) wavelength of this fluorescence transition
        lambda = scatinfo->lambda;
        lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfkobs, scatinfo->velocity);
        return true;
    }

    // resonant scattering
    else if (scatinfo->species < _numElec + _numFluo + _numLine)
    {
        if (scatinfo->dipole)
            _dpf->peeloffScattering(I, Q, U, V, pp->direction(), bfkobs, bfky, pp);
        else
            // isotropic scattering removes polarization, so the contribution is trivially 1
            I = 1.;

        // if coherent (no branching)
        if (scatinfo->lambda == 0.)
            lambda = LyUtils::shiftWavelength(lambda, scatinfo->velocity, pp->direction(), bfkobs);
        // if incoherent (branching)
        else
            lambda = PhotonPacket::shiftedEmissionWavelength(scatinfo->lambda, bfkobs, scatinfo->velocity);

        return false;
    }

    // line emission following FAC inner-shell ionization
    else
    {
        I = 1.;
        lambda = PhotonPacket::shiftedEmissionWavelength(scatinfo->lambda, bfkobs, scatinfo->velocity);
        return true;
    }
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMix::performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const
{
    // draw a random fluorescence channel and atom velocity, unless a previous peel-off stored this already
    setScatteringInfoIfNeeded(const_cast<PhotonPacket*>(pp), state, lambda);
    auto scatinfo = const_cast<PhotonPacket*>(pp)->getScatteringInfo();

    // room for the outgoing direction
    Direction bfknew;

    // Compton scattering, with support for polarization if enabled:
    // determine the new propagation direction and wavelength, and if polarized, update the stokes vector
    if (scatinfo->species < _numElec)
    {
        int i = scatinfo->species;
        const auto& ion = _ionParamv[i];
        lambda = PhotonPacket::shiftedReceptionWavelength(lambda, pp->direction(), scatinfo->velocity);
        bfknew = _com->performScattering(lambda, ion.Z, ion.N, pp->direction(), pp);
        lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfknew, scatinfo->velocity);
    }

    // fluorescence, always unpolarized and isotropic
    else if (scatinfo->species < _numElec + _numFluo)
    {
        // clear the stokes vector (only relevant if polarization support is enabled)
        pp->setUnpolarized();

        // update the photon packet wavelength to the (possibly sampled) wavelength of this fluorescence transition
        lambda = scatinfo->lambda;
        bfknew = random()->direction();
        lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfknew, scatinfo->velocity);

        // indicate that this packet emulates secondary emission;
        pp->setEmulatedSecondaryOrigin(state->mediumIndex());
    }

    // resonant scattering
    else if (scatinfo->species < _numElec + _numFluo + _numLine)
    {
        if (scatinfo->dipole)
        {
            bfknew = _dpf->performScattering(pp->direction(), pp);
        }
        else
        {
            bfknew = random()->direction();
            pp->setUnpolarized();
        }

        // if coherent (no branching)
        if (scatinfo->lambda == 0.)
            lambda = LyUtils::shiftWavelength(lambda, scatinfo->velocity, pp->direction(), bfknew);
        // if incoherent (branching)
        else
            lambda = PhotonPacket::shiftedEmissionWavelength(scatinfo->lambda, bfknew, scatinfo->velocity);
    }

    // line emission following FAC inner-shell ionization, always unpolarized and isotropic
    else
    {
        pp->setUnpolarized();
        lambda = scatinfo->lambda;
        bfknew = random()->direction();
        lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfknew, scatinfo->velocity);
        pp->setEmulatedSecondaryOrigin(state->mediumIndex());
    }

    // execute the scattering event in the photon packet
    pp->scatter(bfknew, state->bulkVelocity(), lambda);
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::indicativeTemperature(const MaterialState* /*state*/, const Array& /*Jv*/) const
{
    return temperature();
}

////////////////////////////////////////////////////////////////////
