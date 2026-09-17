/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "XRayAtomicGasMix.hpp"
#include "ComptonPhaseFunction.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "DipolePhaseFunction.hpp"
#include "FatalError.hpp"
#include "MaterialState.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "Range.hpp"
#include "TextInFile.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // ---- common helper functions ----

    // convert photon energy in eV to and from wavelength in m (same conversion in both directions)
    constexpr double wavelengthToFromEnergy(double x)
    {
        constexpr double front = Constants::h() * Constants::c() / Constants::Qelectron();
        return front / x;
    }

    // return thermal velocity for given gas temperature (in K) and particle mass (in amu)
    double vtherm(double T, double amu)
    {
        return sqrt(Constants::k() / Constants::amu() * T / amu);
    }

    // ---- hardcoded configuration constants ----

    // number of supported atoms; the data provided in the resource files must match this number
    constexpr size_t numAtoms = 30;

    // wavelength range over which our cross sections may be nonzero
    constexpr Range nonZeroRange(wavelengthToFromEnergy(500e3), wavelengthToFromEnergy(4.3));

    // number of wavelengths per dex in high-resolution grid
    constexpr size_t numWavelengthsPerDex = 2500;

    // ---- photo-absorption and fluorescence resources ----

    // basic atom information
    struct AtomParams
    {
        AtomParams(const Array& a) : mass(a[0]), abund(a[1]) {}
        double mass;   // atom mass (amu)
        double abund;  // the default or configured relative abundance (1)
    };

    // photo-absorption cross section parameters
    struct CrossSectionParams
    {
        CrossSectionParams(const Array& a)
            : Z(a[0]), n(a[1]), l(a[2]), Eth(a[3]), Emax(a[4]), E0(a[5]), sigma0(a[6]), ya(a[7]), P(a[8]), yw(a[9]),
              y0(a[10]), y1(a[11])
        {}
        short Z;        // atomic number
        short n;        // principal quantum number of the shell
        short l;        // orbital quantum number of the subshell
        double Eth;     // subshell ionization threshold energy (eV)
        double Emax;    // maximum energy for validity of the formula (eV)
        double E0;      // fit parameter (eV)
        double sigma0;  // fit parameter (Mb = 10^-22 m^2)
        double ya;      // fit parameter (1)
        double P;       // fit parameter (1)
        double yw;      // fit parameter (1)
        double y0;      // fit parameter (1)
        double y1;      // fit parameter (1)
    };

    // fluorescence parameters
    struct FluorescenceParams
    {
        FluorescenceParams(const Array& a) : Z(a[0]), n(a[1]), l(a[2]), omega(a[3]), E(a[4]), W(a[5]) {}
        short Z;       // atomic number
        short n;       // principal quantum number of the shell
        short l;       // orbital quantum number of the subshell
        double omega;  // fluorescence yield (1)
        double E;      // (central) energy of the emitted photon (eV)
        double W;      // FWHM of the Lorentz shape for the emitted photon (eV), or zero
    };

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

    // ---- photo-absorption cross section ----

    // return photo-absorption cross section in m2 for given energy in eV and cross section parameters,
    // without taking into account thermal dispersion
    double crossSection(double E, const CrossSectionParams& p)
    {
        if (E < p.Eth || E >= p.Emax) return 0.;

        double x = E / p.E0 - p.y0;
        double y = std::sqrt(x * x + p.y1 * p.y1);
        double xm1 = x - 1.;
        double Q = 5.5 + p.l - 0.5 * p.P;
        double F = (xm1 * xm1 + p.yw * p.yw) * std::pow(y, -Q) * std::pow(1. + std::sqrt(y / p.ya), -p.P);
        return 1e-22 * p.sigma0 * F;  // from Mb to m2
    }

    // return photo-absorption cross section in m2 for given energy in eV and cross section parameters,
    // approximating thermal dispersion by replacing the steep threshold transition by a sigmoid
    // error function with given parameters (dispersion and maximum value)
    double crossSection(double E, std::pair<double, double> sigmoid, const CrossSectionParams& p)
    {
        double Es, sigmamax;
        std::tie(Es, sigmamax) = sigmoid;
        if (E <= p.Eth - 2. * Es) return 0.;
        if (E >= p.Eth + 2. * Es) return crossSection(E, p);
        return sigmamax * (0.5 + 0.5 * std::erf((E - p.Eth) / Es));
    }
}

////////////////////////////////////////////////////////////////////

void XRayAtomicGasMix::setupSelfBefore()
{
    MaterialMix::setupSelfBefore();
    auto config = find<Configuration>();

    // ---- load resources ----

    // load the atom masses and default abundancies
    auto atomv = loadStruct<AtomParams, 2>(this, "XRay_MA.txt", "atom masses");

    // if the configured abundancies list is nonempty, use it instead of the defaults
    if (!_abundancies.empty())
    {
        if (_abundancies.size() != numAtoms)
            throw FATALERROR("The abundancies list must have exactly " + std::to_string(numAtoms) + " values");
        for (size_t i = 0; i != numAtoms; ++i) atomv[i].abund = _abundancies[i];
    }

    // load the photo-absorption cross section parameters
    auto crossSectionParams = loadStruct<CrossSectionParams, 12>(this, "XRay_PA.txt", "photo-absorption data");

    // load the fluorescence parameters
    auto fluorescenceParams = loadStruct<FluorescenceParams, 6>(this, "XRay_FL.txt", "fluorescence data");

    // create scattering helpers depending on the user-configured implementation type;
    // the respective helper constructors load the required bound-electron scattering resources
    switch (scatterBoundElectrons())
    {
        case BoundElectrons::None:
            _ray = new NoScatteringHelper(this);
            _com = new NoScatteringHelper(this);
            break;
        case BoundElectrons::Free:
            _ray = new NoScatteringHelper(this);
            _com = new FreeComptonHelper(this);
            break;
        case BoundElectrons::FreeWithPolarization:
            _ray = new NoScatteringHelper(this);
            _com = new FreeComptonWithPolarizationHelper(this);
            break;
        case BoundElectrons::Good:
            _ray = new SmoothRayleighHelper(this);
            _com = new BoundComptonHelper(this);
            break;
        case BoundElectrons::Exact:
            _ray = new AnomalousRayleighHelper(this);
            _com = new BoundComptonHelper(this);
            break;
    }

    // ---- thermal dispersion ----

    // calculate the parameters for the sigmoid function approximating the convolution with a Gaussian
    // at the threshold energy for each cross section record, and store the result into a temporary vector;
    // the information includes the thermal energy dispersion at the threshold energy and
    // the intrinsic cross section at the threshold energy plus twice this energy dispersion
    vector<std::pair<double, double>> sigmoidv;
    sigmoidv.reserve(crossSectionParams.size());
    for (const auto& params : crossSectionParams)
    {
        double Es = params.Eth * vtherm(temperature(), atomv[params.Z - 1].mass) / Constants::c();
        double sigmamax = crossSection(params.Eth + 2. * Es, params);
        sigmoidv.emplace_back(Es, sigmamax);
    }

    // calculate and store the thermal velocities corresponding to the scattering channels
    // (Rayleigh scattering for each atom, Compton scattering for each atom, fluorescence for each transition)
    _vthermscav.reserve(2 * numAtoms + fluorescenceParams.size());
    for (int i = 0; i != 2; ++i)
        for (const auto& atom : atomv) _vthermscav.push_back(vtherm(temperature(), atom.mass));
    for (const auto& params : fluorescenceParams)
        _vthermscav.push_back(vtherm(temperature(), atomv[params.Z - 1].mass));

    // ---- wavelength grid ----

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
    int index = 0;
    for (const auto& params : crossSectionParams)
    {
        double Es = sigmoidv[index++].first;
        for (double delta : {-2., -4. / 3., -2. / 3., 0., 2. / 3., 4. / 3., 2.})
        {
            double lambda = wavelengthToFromEnergy(params.Eth + delta * Es);
            if (range.contains(lambda)) lambdav.push_back(lambda);
        }
    }

    // add the fluorescence emission wavelengths
    for (const auto& params : fluorescenceParams)
    {
        double lambda = wavelengthToFromEnergy(params.E);
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

    // ---- extinction ----

    // calculate the extinction cross section at every wavelength; to guarantee that the cross section is zero
    // for wavelengths outside our range, leave the values for the three outer wavelength points at zero
    _sigmaextv.resize(numLambda);
    for (int ell = 1; ell < numLambda - 2; ++ell)
    {
        double lambda = lambdav[ell];
        double sigma = 0.;

        // bound electron scattering
        for (size_t Z = 1; Z <= numAtoms; ++Z)
        {
            sigma += (_ray->sectionSca(lambda, Z, Z) + _com->sectionSca(lambda, Z, Z)) * atomv[Z - 1].abund;
        }
        // enhancement due to coherent scattering on molecular H2 w.r.t. Rayleigh scattering on atomic H
        sigma += _ray->sectionSca(lambda, 1, 1) * atomv[1 - 1].abund * _H2fraction;

        // photo-absorption and fluorescence
        int index = 0;
        for (const auto& params : crossSectionParams)
        {
            double E = wavelengthToFromEnergy(lambda);
            sigma += crossSection(E, sigmoidv[index++], params) * atomv[params.Z - 1].abund;
        }
        _sigmaextv[ell] = sigma;
    }

    // ---- scattering ----

    // calculate and store the fluorescence emission parameters
    int numFluos = fluorescenceParams.size();
    _lambdafluov.resize(numFluos);
    _centralfluov.resize(numFluos);
    _widthfluov.resize(numFluos);
    for (int k = 0; k != numFluos; ++k)
    {
        const auto& params = fluorescenceParams[k];
        if (!params.W)
        {
            // if the line shape has zero width, convert the central enery to the fixed wavelength
            _lambdafluov[k] = wavelengthToFromEnergy(params.E);
        }
        else
        {
            // if the line shape has nonzero width, copy the parameters of the line shape
            // so that we can sample from the Lorentz distribution in energy space
            _centralfluov[k] = params.E;
            _widthfluov[k] = params.W / 2.;  // convert from FWHM to HWHM
        }
    }

    // make room for the scattering cross section and the cumulative fluorescence/scattering probabilities
    _sigmascav.resize(numLambda);
    _cumprobscavv.resize(numLambda, 0);

    // provide temporary array for the non-normalized fluorescence/scattering contributions (at the current wavelength)
    Array contribv(2 * numAtoms + fluorescenceParams.size());

    // calculate the above for every wavelength; as before, leave the values for the outer wavelength points at zero
    for (int ell = 1; ell < numLambda - 2; ++ell)
    {
        double lambda = lambdav[ell];
        double E = wavelengthToFromEnergy(lambda);

        // bound electron scattering
        for (size_t Z = 1; Z <= numAtoms; ++Z)
        {
            contribv[Z - 1] = _ray->sectionSca(lambda, Z, Z) * atomv[Z - 1].abund;
            contribv[numAtoms + Z - 1] = _com->sectionSca(lambda, Z, Z) * atomv[Z - 1].abund;
        }
        // enhancement due to coherent scattering on molecular H2 w.r.t. Rayleigh scattering on atomic H
        contribv[1 - 1] += _ray->sectionSca(lambda, 1, 1) * atomv[1 - 1].abund * _H2fraction;

        // fluorescence: iterate over both cross section and fluorescence parameter sets in sync
        auto flp = fluorescenceParams.begin();
        int index = 0;
        for (const auto& csp : crossSectionParams)
        {
            auto sigmoid = sigmoidv[index++];

            // process all fluorescence parameter sets matching this cross section set
            while (flp != fluorescenceParams.end() && flp->Z == csp.Z && flp->n == csp.n && flp->l == csp.l)
            {
                double contribution = crossSection(E, sigmoid, csp) * atomv[csp.Z - 1].abund * flp->omega;
                contribv[2 * numAtoms + flp - fluorescenceParams.begin()] = contribution;
                flp++;
            }
        }

        // determine the normalized cumulative probability distribution and the cross section
        _sigmascav[ell] = NR::cdf(_cumprobscavv[ell], contribv);
    }
}

////////////////////////////////////////////////////////////////////

XRayAtomicGasMix::~XRayAtomicGasMix()
{
    delete _ray;
    delete _com;
}

////////////////////////////////////////////////////////////////////

int XRayAtomicGasMix::indexForLambda(double lambda) const
{
    return NR::locateClip(_lambdav, lambda);
}

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType XRayAtomicGasMix::materialType() const
{
    return MaterialType::Gas;
}

////////////////////////////////////////////////////////////////////

bool XRayAtomicGasMix::hasPolarizedScattering() const
{
    return scatterBoundElectrons() == BoundElectrons::FreeWithPolarization;
}

////////////////////////////////////////////////////////////////////

bool XRayAtomicGasMix::hasScatteringDispersion() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool XRayAtomicGasMix::scatteringEmulatesSecondaryEmission() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> XRayAtomicGasMix::specificStateVariableInfo() const
{
    return vector<StateVariable>{StateVariable::numberDensity()};
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::mass() const
{
    return Constants::Mproton();
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::sectionAbs(double lambda) const
{
    int index = indexForLambda(lambda);
    return _sigmaextv[index] - _sigmascav[index];
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::sectionSca(double lambda) const
{
    return _sigmascav[indexForLambda(lambda)];
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::sectionExt(double lambda) const
{
    return _sigmaextv[indexForLambda(lambda)];
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double number = state->numberDensity();
    return number > 0. ? sectionAbs(lambda) * number : 0.;
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::opacitySca(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double number = state->numberDensity();
    return number > 0. ? sectionSca(lambda) * number : 0.;
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double number = state->numberDensity();
    return number > 0. ? sectionExt(lambda) * number : 0.;
}

////////////////////////////////////////////////////////////////////

void XRayAtomicGasMix::setScatteringInfoIfNeeded(PhotonPacket::ScatteringInfo* scatinfo, double lambda) const
{
    if (!scatinfo->valid)
    {
        scatinfo->valid = true;
        scatinfo->species = NR::locateClip(_cumprobscavv[indexForLambda(lambda)], random()->uniform());
        if (temperature() > 0.) scatinfo->velocity = _vthermscav[scatinfo->species] * random()->maxwell();

        // for a fluorescence transition, determine the outgoing wavelength from the corresponding parameters
        if (scatinfo->species >= static_cast<int>(2 * numAtoms))
        {
            int i = scatinfo->species - 2 * numAtoms;
            if (_lambdafluov[i])
            {
                // for a zero-width line, simply copy the central wavelength
                scatinfo->lambda = _lambdafluov[i];
            }
            else
            {
                // otherwise sample a wavelength from the Lorentz line shape in energy space;
                // the tails of the Lorentz distribition are very long, occasionaly resulting in negative energies;
                // therefore we loop until the sampled wavelength is meaningful
                while (true)
                {
                    scatinfo->lambda = wavelengthToFromEnergy(_centralfluov[i] + _widthfluov[i] * random()->lorentz());
                    if (nonZeroRange.contains(scatinfo->lambda)) break;
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

bool XRayAtomicGasMix::peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs,
                                         Direction bfky, const MaterialState* /*state*/, const PhotonPacket* pp) const
{
    // draw a random scattering channel and atom velocity, unless a previous peel-off stored this already
    auto scatinfo = const_cast<PhotonPacket*>(pp)->getScatteringInfo();
    setScatteringInfoIfNeeded(scatinfo, lambda);

    // if we have dispersion, for electron scattering, adjust the incoming wavelength to the electron rest frame
    if (temperature() > 0. && scatinfo->species < static_cast<int>(2 * numAtoms))
        lambda = PhotonPacket::shiftedReceptionWavelength(lambda, pp->direction(), scatinfo->velocity);

    // Rayleigh scattering in electron rest frame; no support for polarization
    if (scatinfo->species < static_cast<int>(numAtoms))
    {
        int Z = scatinfo->species + 1;
        _ray->peeloffScattering(I, lambda, Z, Z, pp->direction(), bfkobs);
    }

    // Compton scattering in electron rest frame; with support for polarization if enabled
    else if (scatinfo->species < static_cast<int>(2 * numAtoms))
    {
        int Z = scatinfo->species - numAtoms + 1;
        _com->peeloffScattering(I, Q, U, V, lambda, Z, Z, pp->direction(), bfkobs, bfky, pp);
    }

    // fluorescence
    else
    {
        // unpolarized isotropic emission; the bias weight is trivially 1 and there is no contribution to Q, U, V
        I = 1.;

        // update the photon packet wavelength to the (possibly sampled) wavelength of this fluorescence transition
        lambda = scatinfo->lambda;
    }

    // if we have dispersion, Doppler-shift the outgoing wavelength from the electron rest frame
    if (temperature() > 0.) lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfkobs, scatinfo->velocity);

    // return true for fluorescence, false otherwise
    return scatinfo->species >= static_cast<int>(2 * numAtoms);
}

////////////////////////////////////////////////////////////////////

void XRayAtomicGasMix::performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const
{
    // draw a random fluorescence channel and atom velocity, unless a previous peel-off stored this already
    auto scatinfo = pp->getScatteringInfo();
    setScatteringInfoIfNeeded(scatinfo, lambda);

    // if we have dispersion, for electron scattering, adjust the incoming wavelength to the electron rest frame
    if (temperature() > 0. && scatinfo->species < static_cast<int>(2 * numAtoms))
        lambda = PhotonPacket::shiftedReceptionWavelength(lambda, pp->direction(), scatinfo->velocity);

    // room for the outgoing direction
    Direction bfknew;

    // Rayleigh scattering, no support for polarization: determine the new propagation direction
    if (scatinfo->species < static_cast<int>(numAtoms))
    {
        int Z = scatinfo->species + 1;
        bfknew = _ray->performScattering(lambda, Z, Z, pp->direction());
    }

    // Compton scattering, with support for polarization if enabled:
    // determine the new propagation direction and wavelength, and if polarized, update the stokes vector
    else if (scatinfo->species < static_cast<int>(2 * numAtoms))
    {
        int Z = scatinfo->species - numAtoms + 1;
        bfknew = _com->performScattering(lambda, Z, Z, pp->direction(), pp);
    }

    // fluorescence, always unpolarized and isotropic
    else
    {
        // update the photon packet wavelength to the (possibly sampled) wavelength of this fluorescence transition
        lambda = scatinfo->lambda;

        // draw a random, isotropic outgoing direction
        bfknew = random()->direction();

        // clear the stokes vector (only relevant if polarization support is enabled)
        pp->setUnpolarized();

        // indicate that this packet emulates secondary emission;
        pp->setEmulatedSecondaryOrigin(state->mediumIndex());
    }

    // if we have dispersion, Doppler-shift the outgoing wavelength from the electron rest frame
    if (temperature() > 0.) lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfknew, scatinfo->velocity);

    // execute the scattering event in the photon packet
    pp->scatter(bfknew, state->bulkVelocity(), lambda);
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::indicativeTemperature(const MaterialState* /*state*/, const Array& /*Jv*/) const
{
    return temperature();
}

////////////////////////////////////////////////////////////////////
