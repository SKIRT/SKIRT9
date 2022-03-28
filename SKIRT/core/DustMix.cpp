/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustMix.hpp"
#include "Configuration.hpp"
#include "Log.hpp"
#include "MaterialState.hpp"
#include "NR.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

// built-in constants determining the resolution for discretizing the scattering angles
namespace
{
    // theta from 0 to pi, index t
    constexpr int numTheta = 361;
    constexpr int maxTheta = numTheta - 1;
    constexpr double deltaTheta = M_PI / maxTheta;

    // phi from 0 to 2 pi, index f
    constexpr int numPhi = 361;
    constexpr int maxPhi = numPhi - 1;
    constexpr double deltaPhi = 2 * M_PI / (maxPhi - 1);
}

// hard-coded parameters configuring extreme forward/backward Henyey-Greenstein scattering
namespace
{
    // calculating the value of and sampling from the HG phase function is numerically unstable for abs(g) > gmax
    constexpr double gmax = 0.999999;

    // peel-off bias factors are averaged for abs(g) > glarge
    constexpr double glarge = 0.95;

    // half-opening angle over which to average the HG phase function for abs(g) > glarge
    constexpr double delta = 4. * M_PI / 180.;  // 4 degrees
}

////////////////////////////////////////////////////////////////////

void DustMix::setupSelfAfter()
{
    MaterialMix::setupSelfAfter();
    auto config = find<Configuration>();

    // construct a wavelength grid for sampling dust properties containing two (merged) sets of grid points:
    //  - a fine grid (in log space) covering the entire wavelength range of the simulation
    //  - all specific wavelengths mentioned in the configuration of the simulation (grids, normalizations, ...)
    //    ensuring that the dust properties are interpolated at exactly these wavelengths (e.g. for benchmarks)
    // we first gather all the wavelength points, in arbitrary order, and then sort them
    vector<double> wavelengths;
    wavelengths.reserve(10000);

    // add a fine grid covering the wavelength range of the simulation in log space;
    // use integer multiples as logarithmic grid points so that the grid is stable for changing wavelength ranges
    const double numWavelengthsPerDex = 1000;  // declared as double to avoid casting, but should be an integer
    Range range = config->simulationWavelengthRange();
    int minLambdaSerial = std::floor(numWavelengthsPerDex * log10(range.min()));
    int maxLambdaSerial = std::ceil(numWavelengthsPerDex * log10(range.max()));
    for (int k = minLambdaSerial; k <= maxLambdaSerial; ++k) wavelengths.push_back(pow(10., k / numWavelengthsPerDex));

    // add the wavelengths mentioned in the configuration of the simulation
    for (double lambda : config->simulationWavelengths()) wavelengths.push_back(lambda);

    // sort the wavelengths and remove duplicates
    NR::unique(wavelengths);

    // remove wavelengths beyond 10 cm, if any, because those dust cross sections should be assumed to be zero
    const double dm = 0.1;
    bool radioCutoff = wavelengths.back() > dm;
    if (radioCutoff)
    {
        wavelengths.resize(NR::locate(wavelengths, dm) + 1);
        if (wavelengths.empty() || wavelengths.back() != dm) wavelengths.push_back(dm);
        wavelengths.push_back(dm * 1.001);  // add a dummy border value which is never used
    }

    // remember the wavelength range
    _required.set(wavelengths.front(), wavelengths.back());

    // copy the wavelengths into a temporary (local) array that will be used to sample the dust properties
    Array lambdav = NR::array(wavelengths);
    int numLambda = lambdav.size();

    // derive a wavelength grid that will be used for converting a wavelength to an index in the above array;
    // the grid points are shifted to the left of the actual sample points to approximate rounding
    _lambdav.resize(numLambda);
    _lambdav[0] = lambdav[0];
    for (int ell = 1; ell != numLambda; ++ell)
    {
        _lambdav[ell] = sqrt(lambdav[ell] * lambdav[ell - 1]);
    }

    // get the scattering mode advertised by this dust mix
    auto mode = scatteringMode();

    // if needed, build a scattering angle grid
    if (mode == ScatteringMode::MaterialPhaseFunction || mode == ScatteringMode::SphericalPolarization
        || mode == ScatteringMode::SpheroidalPolarization)
    {
        _thetav.resize(numTheta);
        for (int t = 0; t != numTheta; ++t) _thetav[t] = t * deltaTheta;
    }

    // resize the optical property arrays and tables as needed
    _sigmaabsv.resize(numLambda);
    _sigmascav.resize(numLambda);
    _sigmaextv.resize(numLambda);
    _asymmparv.resize(numLambda);
    if (mode == ScatteringMode::MaterialPhaseFunction || mode == ScatteringMode::SphericalPolarization
        || mode == ScatteringMode::SpheroidalPolarization)
    {
        _S11vv.resize(numLambda, numTheta);
        if (mode == ScatteringMode::SphericalPolarization || mode == ScatteringMode::SpheroidalPolarization)
        {
            _S12vv.resize(numLambda, numTheta);
            _S33vv.resize(numLambda, numTheta);
            _S34vv.resize(numLambda, numTheta);
        }
        if (mode == ScatteringMode::SpheroidalPolarization)
        {
            _sigmaabsvv.resize(numLambda, numTheta);
            _sigmaabspolvv.resize(numLambda, numTheta);
        }
    }

    // obtain the optical properties from the subclass
    _mu = getOpticalProperties(lambdav, _thetav, _sigmaabsv, _sigmascav, _asymmparv, _S11vv, _S12vv, _S33vv, _S34vv,
                               _sigmaabsvv, _sigmaabspolvv);

    // ensure that values for the asymmetry parameter g are not too close to 1 or -1
    int lastAdjustedEll = -1;  // index of longest wavelength for which adjustment has been made
    for (int ell = 0; ell != numLambda; ++ell)
    {
        if (abs(_asymmparv[ell]) > gmax)
        {
            _asymmparv[ell] = std::copysign(gmax, _asymmparv[ell]);
            lastAdjustedEll = ell;
        }
    }
    if (lastAdjustedEll >= 0)
    {
        auto units = find<Units>();
        auto log = find<Log>();
        log->warning(type() + " extreme forward/backward scattering asymmetry parameter values reduced to "
                     + StringUtils::toString(gmax));
        log->warning("  For some or all wavelengths short of "
                     + StringUtils::toString(units->owavelength(lambdav[lastAdjustedEll]), 'g', 3) + " "
                     + units->uwavelength());
    }

    // calculate derived basic optical properties
    for (int ell = 0; ell != numLambda; ++ell)
    {
        _sigmaextv[ell] = _sigmaabsv[ell] + _sigmascav[ell];
    }

    // precalculate discretizations related to the scattering angles as needed
    if (mode == ScatteringMode::MaterialPhaseFunction || mode == ScatteringMode::SphericalPolarization
        || mode == ScatteringMode::SpheroidalPolarization)
    {
        // create a table with the normalized cumulative distribution of theta for each wavelength
        _thetaXvv.resize(numLambda, 0);
        for (int ell = 0; ell != numLambda; ++ell)
        {
            NR::cdf(_thetaXvv[ell], maxTheta, [this, ell](int t) { return _S11vv(ell, t + 1) * sin(_thetav[t + 1]); });
        }

        // create a table with the phase function normalization factor for each wavelength
        _pfnormv.resize(numLambda);
        for (int ell = 0; ell != numLambda; ++ell)
        {
            double sum = 0.;
            for (int t = 0; t != numTheta; ++t)
            {
                sum += _S11vv(ell, t) * sin(_thetav[t]) * deltaTheta;
            }
            _pfnormv[ell] = 2.0 / sum;
        }

        // create tables listing phi, phi/(2 pi), sin(2 phi) and 1-cos(2 phi) for each phi index
        if (mode == ScatteringMode::SphericalPolarization || mode == ScatteringMode::SpheroidalPolarization)
        {
            _phiv.resize(numPhi);
            _phi1v.resize(numPhi);
            _phisv.resize(numPhi);
            _phicv.resize(numPhi);
            for (int f = 0; f != numPhi; ++f)
            {
                double phi = f * deltaPhi;
                _phiv[f] = phi;
                _phi1v[f] = phi / (2 * M_PI);
                _phisv[f] = sin(2 * phi);
                _phicv[f] = 1 - cos(2 * phi);
            }
        }
    }

    // precalculate information to accelerate solving the energy balance equation for the temperature;
    // this is relevant only if the simulation tracks the radiation field
    if (config->hasPanRadiationField())
    {
        _calc.precalculate(this, lambdav, _sigmaabsv);
    }

    // if the wavelength range was cut off, suppress all cross sections beyond the cutoff wavelength
    // by setting the last two values to zero (the penultimate item is for the cutoff wavelength)
    if (radioCutoff)
    {
        _sigmaabsv[numLambda - 2] = _sigmaabsv[numLambda - 1] = 0.;
        _sigmascav[numLambda - 2] = _sigmascav[numLambda - 1] = 0.;
        _sigmaextv[numLambda - 2] = _sigmaextv[numLambda - 1] = 0.;
    }

    // give the subclass a chance to obtain additional precalculated information
    size_t allocatedBytes = initializeExtraProperties(lambdav);

    // calculate and log allocated memory size
    size_t allocatedSize = 0;
    allocatedSize += _thetav.size();
    allocatedSize += _sigmaabsv.size();
    allocatedSize += _sigmascav.size();
    allocatedSize += _sigmaextv.size();
    allocatedSize += _asymmparv.size();
    allocatedSize += _S11vv.size();
    allocatedSize += _S12vv.size();
    allocatedSize += _S33vv.size();
    allocatedSize += _S34vv.size();
    allocatedSize += _thetaXvv.size();
    allocatedSize += _pfnormv.size();
    allocatedSize += _phiv.size();
    allocatedSize += _phi1v.size();
    allocatedSize += _phisv.size();
    allocatedSize += _phicv.size();
    allocatedSize += _sigmaabsvv.size();
    allocatedSize += _sigmaabspolvv.size();

    allocatedBytes += allocatedSize * sizeof(double) + _calc.allocatedBytes();
    find<Log>()->info(type() + " allocated " + StringUtils::toMemSizeString(allocatedBytes) + " of memory");
}

////////////////////////////////////////////////////////////////////

size_t DustMix::initializeExtraProperties(const Array& /*lambdav*/)
{
    return 0;
}

////////////////////////////////////////////////////////////////////

void DustMix::informAvailableWavelengthRange(Range available)
{
    const double fuzzy = 0.011;  // slightly more than 1% because Configuration class extends the range by 1%
    if (!available.containsFuzzy(_required.min(), fuzzy) || !available.containsFuzzy(_required.max(), fuzzy))
    {
        auto units = find<Units>();
        auto outstring = [units](double wavelength) {
            return StringUtils::toString(units->owavelength(wavelength), 'g', 3) + " " + units->uwavelength();
        };

        auto log = find<Log>();
        log->warning(type() + " dust properties are not available for full simulation wavelength range");
        log->warning("  Required: " + outstring(_required.min()) + " - " + outstring(_required.max()));
        log->warning("  Available: " + outstring(available.min()) + " - " + outstring(available.max()));
    }
}

////////////////////////////////////////////////////////////////////

int DustMix::indexForLambda(double lambda) const
{
    return NR::locateClip(_lambdav, lambda);
}

////////////////////////////////////////////////////////////////////

int DustMix::indexForTheta(double theta) const
{
    int t = 0.5 + theta / deltaTheta;
    return max(0, min(t, maxTheta));
}

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType DustMix::materialType() const
{
    return MaterialType::Dust;
}

////////////////////////////////////////////////////////////////////

bool DustMix::hasContinuumEmission() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> DustMix::specificStateVariableInfo() const
{
    return vector<StateVariable>{StateVariable::numberDensity()};
}

////////////////////////////////////////////////////////////////////

double DustMix::mass() const
{
    return _mu;
}

////////////////////////////////////////////////////////////////////

double DustMix::sectionAbs(double lambda) const
{
    return _sigmaabsv[indexForLambda(lambda)];
}

////////////////////////////////////////////////////////////////////

double DustMix::sectionSca(double lambda) const
{
    return _sigmascav[indexForLambda(lambda)];
}

////////////////////////////////////////////////////////////////////

double DustMix::sectionExt(double lambda) const
{
    return _sigmaextv[indexForLambda(lambda)];
}

////////////////////////////////////////////////////////////////////

double DustMix::asymmpar(double lambda) const
{
    return _asymmparv[indexForLambda(lambda)];
}

////////////////////////////////////////////////////////////////////

double DustMix::opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double n = state->numberDensity();
    return n > 0. ? n * _sigmaabsv[indexForLambda(lambda)] : 0.;
}

////////////////////////////////////////////////////////////////////

double DustMix::opacitySca(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double n = state->numberDensity();
    return n > 0. ? n * _sigmascav[indexForLambda(lambda)] : 0.;
}

////////////////////////////////////////////////////////////////////

double DustMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double n = state->numberDensity();
    return n > 0. ? n * _sigmaextv[indexForLambda(lambda)] : 0.;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // This helper function returns the angle phi between the previous and current scattering planes
    // given the normal to the previous scattering plane and the current and new propagation directions
    // of the photon packet. The function returns a zero angle if the light is unpolarized or when the
    // current scattering event is completely forward or backward.
    double angleBetweenScatteringPlanes(Direction np, Direction kc, Direction kn)
    {
        Vec nc = Vec::cross(kc, kn);
        nc /= nc.norm();
        double cosphi = Vec::dot(np, nc);
        double sinphi = Vec::dot(Vec::cross(np, nc), kc);
        double phi = atan2(sinphi, cosphi);
        if (std::isfinite(phi)) return phi;
        return 0.;
    }
}

////////////////////////////////////////////////////////////////////

namespace
{
    // returns the value of the Henyey-Greenstein phase function;
    // numerically stable for abs(g) <= 1-1e-6
    double valueHG(double g, double costheta)
    {
        double t = 1. + g * g - 2. * g * costheta;
        return (1. - g) * (1. + g) / sqrt(t * t * t);
    }

    // returns the value of the definite integral of the Henyey-Greenstein phase function over the interval
    // between two angles 0 <= alpha < beta <= pi specified as cosines -1 <= cosbeta < cosalpha <= 1;
    // numerically stable for values of 0 << abs(g) <= 1-1e-6
    double integralHG(double g, double cosalpha, double cosbeta)
    {
        double ta = sqrt(1. + g * g - 2. * g * cosalpha);
        double tb = sqrt(1. + g * g - 2. * g * cosbeta);
        double f1 = (1. - g) * (1. + g) / g;
        double f2 = (tb - ta) / (tb * ta);
        return f1 * f2;
    }

    // returns the Henyey-Greenstein phase function value averaged over a hard-coded half-opening angle;
    // numerically stable for values of 0 << abs(g) <= 1-1e-6
    double meanHG(double g, double costheta)
    {
        double theta = acos(costheta);
        double cosalpha = cos(theta - delta);
        double cosbeta = cos(theta + delta);
        if (theta < delta)
            return (integralHG(g, 1., cosalpha) + integralHG(g, 1., cosbeta)) / (2. - cosalpha - cosbeta);
        if (theta > M_PI - delta)
            return (integralHG(g, cosalpha, -1.) + integralHG(g, cosbeta, -1.)) / (2. + cosalpha + cosbeta);
        return integralHG(g, cosalpha, cosbeta) / (cosalpha - cosbeta);
    }
}

////////////////////////////////////////////////////////////////////

void DustMix::peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs,
                                Direction bfky, const MaterialState* /*state*/, const PhotonPacket* pp) const
{
    switch (scatteringMode())
    {
        case DustMix::ScatteringMode::HenyeyGreenstein:
        {
            // calculate the value of the Henyey-Greenstein phase function
            double costheta = Vec::dot(pp->direction(), bfkobs);
            double g = asymmpar(lambda);
            double value = abs(g) > glarge ? meanHG(g, costheta) : valueHG(g, costheta);

            // accumulate the weighted sum in the intensity (no support for polarization in this case)
            I += value;
            break;
        }
        case DustMix::ScatteringMode::MaterialPhaseFunction:
        {
            // calculate the value of the material-specific phase function
            double costheta = Vec::dot(pp->direction(), bfkobs);
            double value = phaseFunctionValueForCosine(lambda, costheta);

            // accumulate the weighted sum in the intensity (no support for polarization in this case)
            I += value;
            break;
        }
        case DustMix::ScatteringMode::SphericalPolarization:
        case DustMix::ScatteringMode::SpheroidalPolarization:
        {
            // calculate the value of the material-specific phase function
            double theta = acos(Vec::dot(pp->direction(), bfkobs));
            double phi = angleBetweenScatteringPlanes(pp->normal(), pp->direction(), bfkobs);
            double value = phaseFunctionValue(lambda, theta, phi, pp);

            // copy the polarization state so we can change it without affecting the incoming photon packet
            StokesVector sv = *pp;

            // rotate the Stokes vector reference direction into the scattering plane
            sv.rotateIntoPlane(pp->direction(), bfkobs);

            // apply the Mueller matrix
            applyMueller(lambda, theta, &sv);

            // rotate the Stokes vector reference direction parallel to the instrument frame y-axis
            // it is given bfkobs because the photon is at this point aimed towards the observer
            sv.rotateIntoPlane(bfkobs, bfky);

            // acumulate the weighted sum of all Stokes components to support polarization
            I += value * sv.stokesI();
            Q += value * sv.stokesQ();
            U += value * sv.stokesU();
            V += value * sv.stokesV();
            break;
        }
    }
}

////////////////////////////////////////////////////////////////////

void DustMix::performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const
{
    // determine the new propagation direction and, if required, update the polarization state of the photon packet
    Direction bfknew;
    switch (scatteringMode())
    {
        case DustMix::ScatteringMode::HenyeyGreenstein:
        {
            // sample a scattering angle from the Henyey-Greenstein phase function
            // handle isotropic scattering separately because the HG sampling procedure breaks down in this case
            double g = asymmpar(lambda);
            if (fabs(g) < 1e-6)
            {
                bfknew = random()->direction();
            }
            else
            {
                double f = ((1.0 - g) * (1.0 + g)) / (1.0 - g + 2.0 * g * random()->uniform());
                double costheta = (1.0 + g * g - f * f) / (2.0 * g);
                bfknew = random()->direction(pp->direction(), costheta);
            }
            break;
        }
        case DustMix::ScatteringMode::MaterialPhaseFunction:
        {
            // sample a scattering angle from the material-specific phase function
            double costheta = generateCosineFromPhaseFunction(lambda);
            bfknew = random()->direction(pp->direction(), costheta);
            break;
        }
        case DustMix::ScatteringMode::SphericalPolarization:
        case DustMix::ScatteringMode::SpheroidalPolarization:
        {
            // sample the angles between the previous and new direction from the material-specific phase function,
            // given the incoming polarization state
            double theta, phi;
            std::tie(theta, phi) = generateAnglesFromPhaseFunction(lambda, pp);

            // rotate the Stokes vector (and the scattering plane) of the photon packet
            pp->rotateStokes(phi, pp->direction());

            // apply Mueller matrix to the Stokes vector of the photon packet
            applyMueller(lambda, theta, pp);

            // rotate the propagation direction in the scattering plane
            Vec newdir = pp->direction() * cos(theta) + Vec::cross(pp->normal(), pp->direction()) * sin(theta);

            // normalize the new direction to prevent degradation
            bfknew = Direction(newdir / newdir.norm());
            break;
        }
    }

    // execute the scattering event in the photon packet
    pp->scatter(bfknew, state->bulkVelocity(), lambda);
}

////////////////////////////////////////////////////////////////////

DustMix::ScatteringMode DustMix::scatteringMode() const
{
    return DustMix::ScatteringMode::HenyeyGreenstein;
}

////////////////////////////////////////////////////////////////////

double DustMix::phaseFunctionValueForCosine(double lambda, double costheta) const
{
    int ell = indexForLambda(lambda);
    int t = indexForTheta(acos(costheta));
    return _pfnormv[ell] * _S11vv(ell, t);
}

////////////////////////////////////////////////////////////////////

double DustMix::generateCosineFromPhaseFunction(double lambda) const
{
    return cos(random()->cdfLinLin(_thetav, _thetaXvv[indexForLambda(lambda)]));
}

////////////////////////////////////////////////////////////////////

double DustMix::phaseFunctionValue(double lambda, double theta, double phi, const StokesVector* sv) const
{
    int ell = indexForLambda(lambda);
    int t = indexForTheta(theta);
    double polDegree = sv->linearPolarizationDegree();
    double polAngle = sv->polarizationAngle();
    return _pfnormv[ell] * (_S11vv(ell, t) + polDegree * _S12vv(ell, t) * cos(2. * (phi - polAngle)));
}

////////////////////////////////////////////////////////////////////

std::pair<double, double> DustMix::generateAnglesFromPhaseFunction(double lambda, const StokesVector* sv) const
{
    int ell = indexForLambda(lambda);

    // sample from the normalized cumulative distribution of theta for this wavelength
    double theta = random()->cdfLinLin(_thetav, _thetaXvv[ell]);
    int t = indexForTheta(theta);

    // construct and sample from the normalized cumulative distribution of phi for this wavelength and theta angle
    double polDegree = sv->linearPolarizationDegree();
    double polAngle = sv->polarizationAngle();
    double PF = polDegree * _S12vv(ell, t) / _S11vv(ell, t) / (4 * M_PI);
    double cos2polAngle = cos(2 * polAngle) * PF;
    double sin2polAngle = sin(2 * polAngle) * PF;
    double phi = random()->cdfLinLin(_phiv, _phi1v + cos2polAngle * _phisv + sin2polAngle * _phicv);

    // return the result
    return std::make_pair(theta, phi);
}

////////////////////////////////////////////////////////////////////

void DustMix::applyMueller(double lambda, double theta, StokesVector* sv) const
{
    int ell = indexForLambda(lambda);
    int t = indexForTheta(theta);
    sv->applyMueller(_S11vv(ell, t), _S12vv(ell, t), _S33vv(ell, t), _S34vv(ell, t));
}

////////////////////////////////////////////////////////////////////

const Array& DustMix::thetaGrid() const
{
    return _thetav;
}

////////////////////////////////////////////////////////////////////

const Array& DustMix::sectionsAbs(double lambda) const
{
    int ell = indexForLambda(lambda);
    return _sigmaabsvv[ell];
}

////////////////////////////////////////////////////////////////////

const Array& DustMix::sectionsAbspol(double lambda) const
{
    int ell = indexForLambda(lambda);
    return _sigmaabspolvv[ell];
}

////////////////////////////////////////////////////////////////////

Array DustMix::emissivity(const Array& Jv) const
{
    return _calc.emissivity(Jv);
}

////////////////////////////////////////////////////////////////////

DisjointWavelengthGrid* DustMix::emissionWavelengthGrid() const
{
    return find<Configuration>()->dustEmissionWLG();
}

////////////////////////////////////////////////////////////////////

Array DustMix::emissionSpectrum(const MaterialState* state, const Array& Jv) const
{
    return state->numberDensity() * emissivity(Jv);
}

////////////////////////////////////////////////////////////////////

double DustMix::indicativeTemperature(const MaterialState* /*state*/, const Array& Jv) const
{
    return Jv.size() ? _calc.equilibriumTemperature(0, Jv) : 0.;
}

////////////////////////////////////////////////////////////////////
