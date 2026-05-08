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
#include <cmath>
#include <tuple>

////////////////////////////////////////////////////////////////////

namespace
{
    // ---- common helper functions ----

    static constexpr int numAtoms = 30;  // maximum atomic number used in this class

    // convert photon energy in eV to and from wavelength in m (same conversion in both directions)
    constexpr double wavelengthToFromEnergy(double x)
    {
        constexpr double front = Constants::h() * Constants::c() / Constants::Qelectron();
        return front / x;
    }

    // ---- hardcoded configuration constants ----

    // wavelength range over which our cross sections may be nonzero
    constexpr Range nonZeroRange(wavelengthToFromEnergy(500e3), wavelengthToFromEnergy(4.3));

    // number of wavelengths per dex in high-resolution grid
    constexpr size_t numWavelengthsPerDex = 2500;

    // load data from resource file with N columns into a vector of structs of type S that can be constructed
    // from an array with N elements, and return that vector
    template<class S, int N> vector<S> loadStruct(const SimulationItem* item, string filename, string description)
    {
        vector<S> result;
        TextInFile infile(item, filename, description, true);
        for (int i = 0; i != N; ++i) infile.addColumn(string());
        Array row;
        while (infile.readRow(row)) result.emplace_back(row);
        return result;
    }

    // resource data for photo-absorption
    struct PhotoAbsorbResource
    {
        PhotoAbsorbResource(const Array& a)
            : Z(a[0]), N(a[1]), n(a[2]), l(a[3]), Eth(a[4]), E0(a[5]), sigma0(a[6]), ya(a[7]), P(a[8]), yw(a[9])
        {}

        int ionIndex{-1};                    // index of the ion
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

    // resource data for Lyman-series
    struct LymanResource
    {
        LymanResource(const Array& a) : Z(a[0]), index(a[1]), lamA(a[2]), lam(a[3]) {}
        int ionIndex{-1};  // index of the ion
        double sprob{-1};  // scatter probability after resonant scattering (<=1)
        double vth{-1};    // sqrt(2) * thermal velocity (m/s)
        int Z;             // atomic number
        int index;         // Lyman index (alpha1/2, alpha3/2, beta1/2, ...)
        double lamA;       // wavelength * Einstein A (m/s)
        double lam;        // wavelength (m)

        double section(double lambda) const
        {
            double a = lamA / (4. * M_PI * vth);

            double g = (index % 2) + 1.;
            return LyUtils::section(lambda, lam, vth, a, g);
        }
    };

    // resource data for Lyman branching (incoherent scattering)
    struct LymanBranchResource
    {
        LymanBranchResource(const Array& a) : Z(a[0]), upper(a[1]), lower(a[2]), prob(a[3]) {}
        int Z;        // atomic number
        int upper;    // upper Lyman index
        int lower;    // lower Lyman index
        double prob;  // branching probability
    };
}

////////////////////////////////////////////////////////////////////

// ---- base class for scattering helpers ----

class XRayIonicGasMix::ScatteringHelper
{
public:
    virtual ~ScatteringHelper() {}

    // return scattering cross section for atom in m2
    virtual double sectionSca(double lambda, int Z) const = 0;

    // peel-off unpolarized scattering event: override this in helpers that don't support polarization
    virtual void peeloffScattering(double& /*I*/, double& /*lambda*/, int /*Z*/, Direction /*bfk*/,
                                   Direction /*bfkobs*/) const
    {
        // default implementation does nothing
    }

    // perform unpolarized scattering event: override this in helpers that don't support polarization
    virtual Direction performScattering(double& /*lambda*/, int /*Z*/, Direction /*bfk*/) const
    {
        // default implementation returns null vector
        return Direction();
    }

    // peel-off polarized scattering event: override this in helpers that do support polarization
    virtual void peeloffScattering(double& I, double& /*Q*/, double& /*U*/, double& /*V*/, double& lambda, int Z,
                                   Direction bfk, Direction bfkobs, Direction /*bfky*/,
                                   const StokesVector* /*sv*/) const
    {
        // default implementation calls unpolarized version
        peeloffScattering(I, lambda, Z, bfk, bfkobs);
    }

    // perform polarized scattering event: override this in helpers that do support polarization
    virtual Direction performScattering(double& lambda, int Z, Direction bfk, StokesVector* /*sv*/) const
    {
        // default implementation calls unpolarized version
        return performScattering(lambda, Z, bfk);
    }
};

////////////////////////////////////////////////////////////////////

// ---- no scattering helper ----

namespace
{
    // this helper does nothing; it is used as a stub in case there is no scattering of a given type
    class NoScatteringHelper : public XRayIonicGasMix::ScatteringHelper
    {
    public:
        NoScatteringHelper(SimulationItem* /*item*/) {}

        double sectionSca(double /*lambda*/, int /*Z*/) const override { return 0.; }
    };
}

////////////////////////////////////////////////////////////////////

// ---- free-electron Compton scattering helper ----

namespace
{
    // transition wavelength from Compton to Thomson scattering
    constexpr double comptonWL = wavelengthToFromEnergy(100.);  // 0.1 keV or 12.4 nm

    // this helper forwards all calls to an external helper class for regular Compton scattering
    // (or Thomson scattering for lower energies, because Compton becomes numerically unstable)
    class FreeComptonHelper : public XRayIonicGasMix::ScatteringHelper
    {
    private:
        ComptonPhaseFunction _cpf;
        DipolePhaseFunction _dpf;

    public:
        FreeComptonHelper(SimulationItem* item)
        {
            auto random = item->find<Random>();
            _cpf.initialize(random);
            _dpf.initialize(random);
        }

        double sectionSca(double lambda, int Z) const override
        {
            double sigma = Z * Constants::sigmaThomson();
            if (lambda < comptonWL) sigma *= _cpf.sectionSca(lambda);
            return sigma;
        }

        void peeloffScattering(double& I, double& lambda, int /*Z*/, Direction bfk, Direction bfkobs) const override
        {
            if (lambda < comptonWL)
            {
                double Q, U, V;
                _cpf.peeloffScattering(I, Q, U, V, lambda, bfk, bfkobs, Direction(), nullptr);
            }
            else
            {
                double Q, U, V;
                _dpf.peeloffScattering(I, Q, U, V, bfk, bfkobs, Direction(), nullptr);
            }
        }

        Direction performScattering(double& lambda, int /*Z*/, Direction bfk) const override
        {
            return lambda < comptonWL ? _cpf.performScattering(lambda, bfk, nullptr)
                                      : _dpf.performScattering(bfk, nullptr);
        }
    };
}

////////////////////////////////////////////////////////////////////

// ---- free-electron Compton with polarization scattering helper ----

namespace
{
    // this helper forwards all calls to an external helper class for Compton scattering
    // (or Thomson scattering for lower energies) with support for polarization
    class FreeComptonWithPolarizationHelper : public XRayIonicGasMix::ScatteringHelper
    {
    private:
        ComptonPhaseFunction _cpf;
        DipolePhaseFunction _dpf;

    public:
        FreeComptonWithPolarizationHelper(SimulationItem* item)
        {
            auto random = item->find<Random>();
            _cpf.initialize(random, true);
            _dpf.initialize(random, true);
        }

        double sectionSca(double lambda, int Z) const override
        {
            double sigma = Z * Constants::sigmaThomson();
            if (lambda < comptonWL) sigma *= _cpf.sectionSca(lambda);
            return sigma;
        }

        void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, int /*Z*/, Direction bfk,
                               Direction bfkobs, Direction bfky, const StokesVector* sv) const override
        {
            lambda < comptonWL ? _cpf.peeloffScattering(I, Q, U, V, lambda, bfk, bfkobs, bfky, sv)
                               : _dpf.peeloffScattering(I, Q, U, V, bfk, bfkobs, bfky, sv);
        }

        Direction performScattering(double& lambda, int /*Z*/, Direction bfk, StokesVector* sv) const override
        {
            return lambda < comptonWL ? _cpf.performScattering(lambda, bfk, sv) : _dpf.performScattering(bfk, sv);
        }
    };
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMix::setupSelfBefore()
{
    MaterialMix::setupSelfBefore();

    auto config = find<Configuration>();

    // ------------ parse user properties ------------

    // parse all required ions from the ions property
    string ionString = StringUtils::squeeze(ions());
    if (ionString.empty()) throw FATALERROR("No ions specified");
    for (string ion : StringUtils::split(ionString, ","))
    {
        int Z, N;
        std::tie(Z, N) = AtomUtils::parseIon(ion);
        _ionParamv.emplace_back(Z, N);
    }
    _numIons = _ionParamv.size();

    if (_numIons != (int)abundances().size()) throw FATALERROR("Number of ions and abundances do not match");

    // create scattering helpers depending on the user-configured implementation type
    switch (electronScattering())
    {
        case ElectronScattering::None: _com = new NoScatteringHelper(this); break;
        case ElectronScattering::Free: _com = new FreeComptonHelper(this); break;
        case ElectronScattering::FreeWithPolarization: _com = new FreeComptonWithPolarizationHelper(this); break;
    }
    if (resonantScattering())
    {
        _dpf = new DipolePhaseFunction();
        _dpf->initialize(random(), true);
    }

    // resources that are maintained during the setup
    vector<PhotoAbsorbResource> usedPar;
    vector<FluorescenceResource> usedFlr;
    vector<LymanResource> usedLyr;

    // Use nested scope to load and preprocess resources and discard unused resources
    {
        // ------------ load full resources ------------

        // photo-absorption data
        auto paResource = loadStruct<PhotoAbsorbResource, 10>(this, "Ionic_PA.txt", "photoabsorption data");
        // fluorescence data
        auto flResource = loadStruct<FluorescenceResource, 7>(this, "Ionic_FL.txt", "fluorescence data");
        // generic Lyman series data
        auto lyResource = loadStruct<LymanResource, 4>(this, "Ionic_LY.txt", "lyman series data");
        // Lyman recombination temperature-dependent yields
        StoredTable<3> lyyResource(this, "Ionic_LY_Y.stab", "Z(1),Ly(1),T(K)", "Y(1)");
        // Lyman branching probabilities
        vector<LymanBranchResource> lybResource;
        if (resonantScattering())
            lybResource = loadStruct<LymanBranchResource, 4>(this, "Ionic_LY_B.txt", "branching probabilities");

        // ------------ preprocess resources ------------

        // Lyman recombination can be modelled as fluorescence following an inner shell PA (n,l)=(1,0)
        // This ignores the cascade and only models the transition back to the inner shell.
        // There is no PA data for (n,l)=(1,0) so we can simply add them without worrying about duplicates.
        flResource.reserve(flResource.size() + lyResource.size());
        for (const auto& lyr : lyResource)
        {
            double Z = lyr.Z;
            double Ly = lyr.index;
            double E = wavelengthToFromEnergy(lyr.lam);
            double omega = lyyResource(Z, Ly, temperature());
            Array params = {Z, 1, 1, 0, omega, E, 0.};  // ensure order is correct here!
            flResource.emplace_back(params);
        }

        // ------------ discard unused resources ------------

        // for each (used) ion save all the photo-absorption, fluorescence and resonant Lyman transitions

        // for each ION
        for (int i = 0; i < _numIons; i++)
        {
            auto& ion = _ionParamv[i];

            // add all PA for this ION
            for (auto& pa : paResource)
            {
                if (pa.Z == ion.Z && pa.N == ion.N)
                {
                    pa.ionIndex = i;
                    usedPar.push_back(pa);

                    // add all FL for this PA
                    for (auto& fl : flResource)
                    {
                        if (fl.Z == pa.Z && fl.N == pa.N && fl.n == pa.n && fl.l == pa.l)
                        {
                            int p = usedPar.size() - 1;
                            fl.paIndex = p;
                            usedFlr.push_back(fl);
                        }
                    }
                }
            }

            // add RS for this (H-like) ION
            if (resonantScattering() && ion.N == 1)
            {
                for (auto& ly : lyResource)
                {
                    if (ly.Z == ion.Z)
                    {
                        ly.ionIndex = i;
                        usedLyr.push_back(ly);
                    }
                }
            }
        }
        _numFluo = usedFlr.size();
        _numLym = usedLyr.size();

        // ------------ postprocess used resources ------------

        // calculate the persistent thermal velocities
        // doesn't actually use any resources, but stores this for all 30 atomic numbers
        _vthermv.resize(numAtoms, 0.);
        for (int Z = 1; Z <= numAtoms; Z++) _vthermv[Z - 1] = sqrt(Constants::k() * temperature() / AtomUtils::mass(Z));

        // Photo-absorption
        // calculate the parameters for the sigmoid function approximating the convolution with a Gaussian
        // at the threshold energy for each cross section record, and store the result into a temporary vector;
        // the information includes the thermal energy dispersion at the threshold energy and
        // the intrinsic cross section at the threshold energy plus twice this energy dispersion
        for (auto& upa : usedPar)
        {
            auto& ion = _ionParamv[upa.ionIndex];
            upa.Es = upa.Eth * vtherm(ion.Z) / Constants::c();
            upa.sigmamax = upa.photoAbsorbSection(upa.Eth + 2. * upa.Es);
        }

        // Resonant scattering
        // Calculate the total branching probability for each resonant Lyman transition.
        // This is the probability that a new photon will be emitted after a 'resonant absorption' event.
        // This value can be lower than 1 because low energy photons are ignored and thus not re-emitted.
        // There is some code duplication since we do this later in the calculating of the persistent data.
        if (resonantScattering())
        {
            for (auto& uly : usedLyr)
            {
                // strore thermal velocity for convenience
                uly.vth = M_SQRT2 * vtherm(uly.Z);

                // total branching probability
                uly.sprob = 0.;
                // add up all the probabilities for each current->lower branches
                for (auto& b : lybResource)
                {
                    if (b.Z == uly.Z && b.upper == uly.index) uly.sprob += b.prob;
                    if (uly.sprob == 1.) break;  // speed up since branching matrix has a lot of 1s and 0s
                }
            }
        }

        // ------------ calculate/store persistent data ------------

        // The persistent data is the data that is needed beyond the setup (scattering)
        // No changes should be made to the usedFlr, usedLyr, or lybr arrays after this point.
        // The usedFlr and usedLyr need to be in the same order as the persistent params!

        // Fluorescence
        // Store the Z, wavelength, and width of each fluorescence transition.
        // These are needed when scattering photons after a photo-absorption event.
        _fluorescenceParamv.resize(_numFluo);
        for (int f = 0; f != _numFluo; ++f)
        {
            const auto& ufl = usedFlr[f];
            auto& flp = _fluorescenceParamv[f];

            flp.Z = ufl.Z;
            flp.lambda = wavelengthToFromEnergy(ufl.E);
            flp.width = ufl.W / 2.;  // convert from FWHM to HWHM
        }

        // Resonant scattering
        // Store the Z, Lyman index, wavelength, Voigt parameter and the cumulative branching for each resonant transition.
        // These are needed to sample atom velocities and to determine the branch to scatter to.
        if (resonantScattering())
        {
            _lymanParamv.resize(_numLym);
            for (int l = 0; l != _numLym; ++l)
            {
                const auto& uly = usedLyr[l];
                auto& lyp = _lymanParamv[l];

                double vth = M_SQRT2 * vtherm(uly.Z);

                lyp.Z = uly.Z;
                lyp.index = uly.index;
                lyp.lambda = uly.lam;
                lyp.a = uly.lamA / (4. * M_PI * vth);

                int Z = lyp.Z;
                int upper = lyp.index;

                // branching probability
                Array pLyl(0., upper + 1);  // can only decay to lower Lyman index
                for (auto& b : lybResource)
                {
                    // for each Lyman Z, upper level, store all the lower probabilities
                    if (Z == b.Z && upper == b.upper) pLyl[b.lower] = b.prob;
                }

                // store the cumulative
                NR::cdf(lyp.cumbranchingv, pLyl);
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

    // we first gather all the wavelength points, in arbitrary order, and then sort them
    vector<double> lambdav;
    lambdav.reserve(5 * numWavelengthsPerDex);

    // get the relevant range (intersection of simulation range and nonzero range)
    Range range = config->simulationWavelengthRange();
    range.intersect(nonZeroRange);

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
    for (const auto& upa : usedPar)
    {
        double Es = upa.Es;
        for (double delta : {-2., -4. / 3., -2. / 3., 0., 2. / 3., 4. / 3., 2.})
        {
            double lambda = wavelengthToFromEnergy(upa.Eth + delta * Es);
            if (range.contains(lambda)) lambdav.push_back(lambda);
        }
    }

    // add the fluorescence emission wavelengths
    for (const auto& ufl : usedFlr)
    {
        double lambda = wavelengthToFromEnergy(ufl.E);
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
        double lambda = lambdav[ell];
        double sigma = 0.;

        // electron scattering
        for (int i = 0; i < _numIons; i++)
        {
            const auto& ion = _ionParamv[i];
            sigma += _com->sectionSca(lambda, ion.Z) * _abundances[i];
        }

        // photo-absorption and fluorescence
        for (const auto& upa : usedPar)
        {
            double E = wavelengthToFromEnergy(lambda);
            sigma += upa.photoAbsorbThermalSection(E) * _abundances[upa.ionIndex];
        }

        // resonant scattering
        for (const auto& uly : usedLyr)
        {
            sigma += uly.section(lambda) * _abundances[uly.ionIndex];
        }

        _sigmaextv[ell] = sigma;
    }

    // ------------ scattering ------------

    // make room for the scattering cross section and the cumulative fluorescence/scattering probabilities
    _sigmascav.resize(numLambda);
    _cumsigmascavv.resize(numLambda, 0);

    // provide temporary array for the non-normalized fluorescence/scattering contributions (at the current wavelength)
    int numInteractions = _numIons + _numFluo + _numLym;
    Array sections(numInteractions);

    // calculate the above for every wavelength; as before, leave the values for the outer wavelength points at zero
    for (int ell = 1; ell < numLambda - 2; ++ell)
    {
        double lambda = lambdav[ell];
        double E = wavelengthToFromEnergy(lambda);

        // electron scattering
        for (int i = 0; i < _numIons; i++)
        {
            const auto& ion = _ionParamv[i];

            sections[i] = _com->sectionSca(lambda, ion.Z) * _abundances[i];
        }

        // fluorescence: iterate over both cross section and fluorescence parameter sets in sync
        for (int f = 0; f < _numFluo; f++)
        {
            const auto& ufl = usedFlr[f];
            const auto& upa = usedPar[ufl.paIndex];

            double section = upa.photoAbsorbThermalSection(E) * _abundances[upa.ionIndex] * ufl.omega;
            sections[_numIons + f] = section;
        }

        // resonant scattering
        for (int l = 0; l < _numLym; l++)
        {
            const auto& uly = usedLyr[l];

            double section = uly.section(lambda) * _abundances[uly.ionIndex] * uly.sprob;
            sections[_numIons + _numFluo + l] = section;
        }

        // determine the normalized cumulative probability distribution and the cross section
        _sigmascav[ell] = NR::cdf(_cumsigmascavv[ell], sections);
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
        if (scatinfo->species < _numIons)
        {
            int i = scatinfo->species;
            const auto& ion = _ionParamv[i];
            scatinfo->velocity = vtherm(ion.Z) * random()->maxwell();
        }
        // Fluorescenct emission (scattering)
        else if (scatinfo->species < _numIons + _numFluo)
        {
            int f = scatinfo->species - _numIons;
            auto flp = _fluorescenceParamv[f];
            double lambda = flp.lambda;
            double width = flp.width;

            scatinfo->velocity = vtherm(flp.Z) * random()->maxwell();

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
        // Resonant Lyman scattering
        else
        {
            int l = scatinfo->species - _numIons - _numFluo;
            const auto& ulyp = _lymanParamv[l];  // upper Lyman

            int upper = ulyp.index;
            int lower = NR::locateFail(ulyp.cumbranchingv, random()->uniform());

            if (lower == -1) throw FATALERROR("Sampling from Lyman branching probability has failed");

            double center;
            double a;
            double vth;
            bool J32;  // true if Ju=3/2 -> happens at odd Lyman index

            // if coherent (no branching)
            if (lower == upper)
            {
                vth = M_SQRT2 * vtherm(ulyp.Z);
                a = ulyp.a;
                center = ulyp.lambda;
                J32 = upper % 2 == 1;

                scatinfo->lambda = 0.;  // explicitly don't use
            }
            // if incoherent (branching)
            else
            {
                int index = l - (upper - lower);  // index of the lower branching
                if (index < 0 || index >= _numLym) throw FATALERROR("upper/lower index out of range");

                const auto& llyp = _lymanParamv[index];  // lower Lyman

                // set parameters to those of the lower branching
                vth = M_SQRT2 * vtherm(llyp.Z);
                a = llyp.a;
                center = llyp.lambda;
                J32 = false;  // branching is isotropic

                scatinfo->lambda = llyp.lambda;
            }

            // if J32 -> Lya1, Lyb1, ... -> 50/50   dipole/isotropic
            // if J12 -> Lya2, Lyb2, ... -> 100     isotropic
            scatinfo->dipole = J32 ? random()->uniform() < 0.5 : false;

            // sample a atom velocity from Voigt profile
            scatinfo->velocity = LyUtils::sampleAtomVelocity(
                lambda, center, vth, a, temperature(), state->numberDensity(), pp->direction(), config(), random());
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
    if (scatinfo->species < _numIons)
    {
        int i = scatinfo->species - _numIons;
        const auto& ion = _ionParamv[i];
        // transform the wavelength into the rest frame of the electron
        lambda = PhotonPacket::shiftedReceptionWavelength(lambda, pp->direction(), scatinfo->velocity);
        _com->peeloffScattering(I, Q, U, V, lambda, ion.Z, pp->direction(), bfkobs, bfky, pp);
        lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfkobs, scatinfo->velocity);
        return false;
    }

    // fluorescence
    else if (scatinfo->species < _numIons + _numFluo)
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
    if (scatinfo->species < _numIons)
    {
        int i = scatinfo->species - _numIons;
        const auto& ion = _ionParamv[i];
        lambda = PhotonPacket::shiftedReceptionWavelength(lambda, pp->direction(), scatinfo->velocity);
        bfknew = _com->performScattering(lambda, ion.Z, pp->direction(), pp);
        lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfknew, scatinfo->velocity);
    }

    // fluorescence, always unpolarized and isotropic
    else if (scatinfo->species < _numIons + _numFluo)
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
