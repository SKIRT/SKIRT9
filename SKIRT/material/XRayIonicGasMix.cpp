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
#include "FatalError.hpp"
#include "LyUtils.hpp"
#include "MaterialState.hpp"
#include "NR.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"
#include "Range.hpp"
#include "StoredTable.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <unordered_set>

////////////////////////////////////////////////////////////////////

namespace
{
    // ---- hardcoded configuration constants ----

    constexpr int numAtoms = 30;         // maximum atomic number used in this class
    constexpr int resourceN[] = {1, 2};  // the electron numbers for the available resources

    // ---- common helper resources----

    // return true if N is in the resourceN list, indicating that the recombination resource for N is present
    constexpr bool hasResource(int N)
    {
        // manual implementation as std::find() is not constexpr in C++14
        for (int r : resourceN)
            if (r == N) return true;
        return false;
    }

    // convert photon energy in eV to and from wavelength in m (same conversion in both directions)
    constexpr double wavelengthToFromEnergy(double x)
    {
        constexpr double front = Constants::h() * Constants::c() / Constants::Qelectron();
        return front / x;
    }

    // calculate the Voigt a parameter
    constexpr double voigtA(double lamGamma, double vth)
    {
        return lamGamma / (4. * M_PI * vth * M_SQRT2);
    }

    // wavelength range over which our cross sections may be nonzero
    constexpr Range nonZeroRange(wavelengthToFromEnergy(500e3), wavelengthToFromEnergy(4.3));

    // number of wavelengths per dex in high-resolution grid
    constexpr size_t numWavelengthsPerDex = 1000;

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

        if (deltaJ > 0)  // 1
        {
            E1 = 0.1 * 3.0 * j * (6.0 * j + 7.0) / ((j + 1.0) * (2.0 * j + 1.0));
            E2 = 0.1 * (2.0 * j + 5.0) * (j + 2.0) / ((j + 1.0) * (2.0 * j + 1.0));
        }
        else if (deltaJ == 0)  // 0
        {
            E1 = 0.1 * 3.0 * (2.0 * j * j + 2.0 * j + 1.0) / (j * (j + 1.0));
            E2 = 0.1 * (2.0 * j - 1.0) * (2.0 * j + 3.0) / (j * (j + 1.0));
        }
        else if (deltaJ < 0)  // -1
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

        double Es;
        double sigmamax;

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
            : Z(a[0]), N(a[1]), lineIndex(a[2]), lowerJ(a[3]), upperJ(a[4]), lam(a[5]), A(a[6]), lamGamma(a[7])
        {}

        int ionIndex{-1};        // index of the corresponding ion in the user-specified ion list
        double scatterProb{-1};  // total probability that resonant absorption leads to re-emission

        int Z;            // atomic number
        int N;            // number of electrons
        int lineIndex;    // line index used to identify this transition
        double lowerJ;    // total angular momentum of the lower level
        double upperJ;    // total angular momentum of the upper level
        double lam;       // transition wavelength (m)
        double A;         // Einstein A coefficient (1/s)
        double lamGamma;  // lambda times natural line width (m)

        double section(double lambda, double vth) const
        {
            double a = voigtA(lamGamma, vth);
            double g = 2.0 * upperJ + 1.0;
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

    // This function loads data from a resource file, only keeping the required data for the present ions.
    // The data is loaded into a vector of structs of type S that can be constructed from an array with C elements,
    // i.e. C columns. The first two columns of the data must always be the ion parameters Z and N.
    // The function uses an set of unique ion hashes to quickly check if an ion is present.
    // These hashes must be calculated from the AtomUtils::ionIndex() function.
    template<class S, int C>
    vector<S> loadPresent(SimulationItem* item, string filename, string description,
                          const std::unordered_set<int>& ionSet)
    {
        vector<S> result;
        TextInFile infile(item, filename, description, true);
        infile.addColumn("Z");
        infile.addColumn("N");
        for (int i = 2; i != C; ++i) infile.addColumn(string());

        Array row;
        while (infile.readRow(row))
        {
            int Z = row[0];
            int N = row[1];
            int hash = AtomUtils::ionIndex(Z, N);

            // only keep ions that are present
            if (ionSet.count(hash)) result.emplace_back(row);
        }
        return result;
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
    /*
	   The setup performs all prior calculations that are needed during the simulation.
	   This is organized in the following steps:
	   1. Parse user properties
	   2. Load the resources for all present ions
	   3. Preprocess resources by adding more fluorescent data that is implied from the resonant data
	   4. Postprocess resources by adding data that is not present in the resources
	   5. Calculate the persistent data that is needed during the simulation.
	   6. Calculate and store the cross section for each wavelength 
	*/

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

    // check if any ions that allow lyman scattering are present
    auto it = std::find_if(usedN.begin(), usedN.end(), hasResource);
    if (resonantScattering() && it == usedN.end())
        throw FATALERROR("Resonant scattering requires ions that allow resonant scattering");

    // create scattering helpers depending on the user-configured implementation type
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

    // photo-absorption data
    auto paResources = loadPresent<PhotoAbsorbResource, 10>(this, "Ionic_PA.txt", "photo-absorption data", ionSet);

    // fluorescence data
    auto fluoResources = loadPresent<FluorescenceResource, 7>(this, "Ionic_FL.txt", "fluorescence data", ionSet);

    // generic line data (recombination and resonant scattering)
    auto lineResources = loadPresent<LineResource, 8>(this, "Ionic_LN.txt", "resonant line data", ionSet);

    // branching probability data
    vector<BranchResource> branchResources;
    if (resonantScattering())
        branchResources =
            loadPresent<BranchResource, 5>(this, "Ionic_BR.txt", "resonant branching probabilities", ionSet);

    // yields for recombination for each value of N
    std::map<int, StoredTable<3>> recoResources;
    for (int N : usedN)
    {
        if (!hasResource(N)) continue;

        // Can't copy StoredTable so must use C++14 emplace
        recoResources.emplace(std::piecewise_construct, std::forward_as_tuple(N),
                              std::forward_as_tuple(this, branchRrFilename(N), "Z(1),Index(1),T(K)", "Y(1)"));
    }

    // ------------ preprocess resources ------------

    // Recombination lines can be modelled as scattering following an inner-shell photo-absorption (n,l)=(1,0).
    // This ignores the cascade and only models the transition back to the inner shell.
    // We can thus model it the same as fluorescence and will add it to the fluorescence resources.
    for (const auto& lineRes : lineResources)
    {
        double Z = lineRes.Z;
        double N = lineRes.N;
        double lineIndex = lineRes.lineIndex;
        double E = wavelengthToFromEnergy(lineRes.lam);
        double omega = recoResources[lineRes.N](Z, lineIndex, temperature());

        if (omega == 0) continue;

        Array params = {Z, N, 1., 0., omega, E, 0.};
        fluoResources.emplace_back(params);
    }

    // Li-like II -> He-like z.
    for (const auto& lineRes : lineResources)
    {
        // Li-like and z lines only
        if (lineRes.N != 3 || lineRes.lineIndex != 0) continue;

        double Z = lineRes.Z;
        double N = lineRes.N;
        double E = wavelengthToFromEnergy(lineRes.lam);
        double omega = 0.75;  // inner-shell ionisation of Li-like populates upper level of z with 3/4

        Array params = {Z, N, 1., 0., omega, E, 0.};
        fluoResources.emplace_back(params);
    }

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
    // No changes should be made to the usedFlr, usedLines, or usedBranchRs arrays after this point.
    // The usedFlr and usedLines arrays need to remain consistent with the persistent parameter arrays.

    // Fluorescence
    // Store the Z, wavelength, and width of each fluorescence transition.
    // These are needed when scattering photons after a photo-absorption event.
    _fluorescenceParamv.resize(_numFluo);
    for (int f = 0; f != _numFluo; ++f)
    {
        const auto& fluoRes = fluoResources[f];
        auto& fluo = _fluorescenceParamv[f];

        fluo.vth = vtherm(fluoRes.Z);
        fluo.lambda = wavelengthToFromEnergy(fluoRes.E);
        fluo.width = fluoRes.W / 2.;  // convert from FWHM to HWHM
    }

    // Resonant scattering
    // Store the ion index, Z, line index, wavelength, Voigt parameter and the cumulative branching
    // for each resonant transition. These are needed to sample atom velocities and to determine
    // the branch to scatter to.
    if (resonantScattering())
    {
        _resonantParamv.resize(_numLine);

        for (int r = 0; r != _numLine; ++r)
        {
            const auto& lineRes = lineResources[r];
            auto& res = _resonantParamv[r];

            res.ionIndex = lineRes.ionIndex;
            res.lineIndex = lineRes.lineIndex;
            res.vth = vtherm(lineRes.Z);
            res.lambda = lineRes.lam;
            res.a = voigtA(lineRes.lamGamma, res.vth);
            res.dipoleFraction = dipoleFractionFromJ(lineRes.lowerJ, lineRes.upperJ);

            int upper = res.lineIndex;

            // branching probability
            Array scatProb(0., upper + 1);  // can only decay to lower index
            for (const auto& braRes : branchResources)
            {
                if (braRes.ionIndex == res.ionIndex && braRes.upperIndex == upper)
                    scatProb[braRes.lowerIndex] = braRes.prob;
            }

            // store the cumulative
            NR::cdf(res.cumBranchingv, scatProb);
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
    for (int k = minLambdaSerial; k <= maxLambdaSerial; ++k) lambdav.push_back(pow(10., k / numPerDex));

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
    for (int ell = 1; ell != numLambda; ++ell)
    {
        _lambdav[ell] = sqrt(lambdav[ell] * lambdav[ell - 1]);
    }

    // ------------ extinction ------------

    // calculate the extinction cross section at every wavelength; to guarantee that the cross section is zero
    // for wavelengths outside our range, leave the values for the three outer wavelength points at zero
    _sigmaextv.resize(numLambda);
    for (int ell = 1; ell < numLambda - 2; ++ell)
    {
        double lambda = _lambdav[ell];
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

        // resonant scattering
        if (resonantScattering())
        {
            for (const auto& lineRes : lineResources)
            {
                double vth = vtherm(lineRes.Z);
                sigma += lineRes.section(lambda, vth) * abundances[lineRes.ionIndex];
            }
        }

        _sigmaextv[ell] = sigma;
    }

    // ------------ scattering ------------

    // make room for the scattering cross section and the cumulative fluorescence/scattering probabilities
    _sigmascav.resize(numLambda);
    _cumsigmascavv.resize(numLambda, 0);

    // provide temporary array for the non-normalized fluorescence/scattering contributions (at the current wavelength)
    int numInteractions = _numElec + _numFluo + _numLine;
    Array sigmas(numInteractions);

    // calculate the above for every wavelength; as before, leave the values for the outer wavelength points at zero
    for (int ell = 1; ell < numLambda - 2; ++ell)
    {
        double lambda = _lambdav[ell];
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
                double sigma = lineRes.section(lambda, vth) * abundances[lineRes.ionIndex] * lineRes.scatterProb;
                sigmas[_numElec + _numFluo + r] = sigma;
            }
        }

        // determine the normalized cumulative probability distribution and the cross section
        _sigmascav[ell] = NR::cdf(_cumsigmascavv[ell], sigmas);
    }
}

////////////////////////////////////////////////////////////////////

XRayIonicGasMix::~XRayIonicGasMix()
{
    delete _com;
    if (resonantScattering()) delete _dpf;
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
        else
        {
            int r = scatinfo->species - _numElec - _numFluo;
            const auto& res = _resonantParamv[r];

            int upper = res.lineIndex;
            int lower = NR::locateFail(res.cumBranchingv, random()->uniform());

            if (lower == -1) throw FATALERROR("Sampling from resonant branching probability has failed");

            double center;
            double a;
            double vth;
            double dipoleFraction;

            // if coherent (no branching)
            if (lower == upper)
            {
                vth = res.vth;
                a = res.a;
                center = res.lambda;
                dipoleFraction = res.dipoleFraction;

                scatinfo->lambda = 0.;  // inform scatter is coherent
            }
            // if incoherent (branching)
            else
            {
                // find  index of the lower resonantParam
                int lr = r - (upper - lower);  // will work since resonantParamv is sorted using Z, N, lineIndex
                if (lr < 0 || lr >= _numLine) throw FATALERROR("upper/lower index out of range");

                const auto& lres = _resonantParamv[lr];

                // // set parameters to those of the lower branching
                vth = lres.vth;
                a = lres.a;
                center = lres.lambda;
                dipoleFraction = 0.0;  // branching is isotropic

                scatinfo->lambda = lres.lambda;
            }

            // determine whether this scattering event uses the dipole phase function
            scatinfo->dipole = random()->uniform() < dipoleFraction;

            // sample an atom velocity from the Voigt profile
            scatinfo->velocity =
                LyUtils::sampleAtomVelocity(lambda, center, M_SQRT2 * vth, a, temperature(), state->numberDensity(),
                                            pp->direction(), config(), random());
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
    else
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
    else
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

    // execute the scattering event in the photon packet
    pp->scatter(bfknew, state->bulkVelocity(), lambda);
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::indicativeTemperature(const MaterialState* /*state*/, const Array& /*Jv*/) const
{
    return temperature();
}

////////////////////////////////////////////////////////////////////
