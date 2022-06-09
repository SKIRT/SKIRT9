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
    double vtherm(double T, double amu) { return sqrt(Constants::k() / Constants::amu() * T / amu); }

    // multiplicator to convert energy in keV to scaled energy E / (m_e c^2)
    constexpr double keVtoScaledEnergy =
        (1e3 * Constants::Qelectron()) / (Constants::Melectron() * Constants::c() * Constants::c());

    // multiplicator to convert scaled energy to energy in units of 12.4 keV
    constexpr double scaledEnergyTo12keV =
        (Constants::Melectron() * Constants::c() * Constants::c()) / (12.4e3 * Constants::Qelectron());

    // convert wavelength to scaled photon energy: h nu / (m_e c^2)
    constexpr double scaledEnergy(double lambda)
    {
        constexpr double front = Constants::h() / Constants::Melectron() / Constants::c();
        return front / lambda;
    }

    // ---- hardcoded configuration constants ----

    // number of supported atoms; the data provided in the resource files must match this number
    constexpr size_t numAtoms = 30;

    // wavelength range over which our cross sections may be nonzero
    constexpr Range nonZeroRange(wavelengthToFromEnergy(500e3), wavelengthToFromEnergy(4.3));

    // number of wavelengths per dex in high-resolution grid
    constexpr size_t numWavelengthsPerDex = 2500;

    // discretization of the phase function over scattering angle: theta from 0 to pi, index t
    constexpr size_t numTheta = 361;
    constexpr size_t maxTheta = numTheta - 1;
    constexpr double deltaTheta = M_PI / maxTheta;

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
        FluorescenceParams(const Array& a) : Z(a[0]), n(a[1]), omega(a[2]), E(a[3]) {}
        short Z;       // atomic number
        short n;       // principal quantum number of the shell (always 1 because we only support K shell)
        double omega;  // fluorescence yield (1)
        double E;      // energy of the emitted photon (eV)
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

    // ---- bound-electron scattering resources ----

    // load data from resource file with N columns into a vector of N arrays, and return that vector;
    // each of the arrays is resized to remove trailing NaN values, if applicable
    vector<Array> loadColumns(int N, const SimulationItem* item, string filename, string description)
    {
        TextInFile infile(item, filename, description, true);
        for (int i = 0; i != N; ++i) infile.addColumn(string());
        vector<Array> columns = infile.readAllColumns();

        // clip any columns with trailing NaNs
        for (Array& column : columns)
        {
            size_t n = column.size();
            while (n && std::isnan(column[n - 1])) --n;
            if (n != column.size())
            {
                // we need to make a copy because resizing an array clears its contents
                Array copy = column;
                column.resize(n);
                for (size_t i = 0; i != n; ++i) column[i] = copy[i];
            }
        }
        return columns;
    }
}

////////////////////////////////////////////////////////////////////

// ---- base class for scattering helpers ----

class XRayAtomicGasMix::ScatteringHelper
{
public:
    virtual ~ScatteringHelper() {}

    // return scattering cross section for atom in m2
    virtual double sectionSca(double lambda, int Z) const = 0;

    // peel-off scattering event
    virtual void peeloffScattering(double& I, double& lambda, int Z, Direction bfk, Direction bfkobs) const = 0;

    // perform scattering event
    virtual Direction performScattering(double& lambda, int Z, Direction bfk) const = 0;
};

////////////////////////////////////////////////////////////////////

// ---- no scattering helper ----

namespace
{
    // this helper does nothing; it is used as a stub in case there is no scattering of a given type
    class NoScatteringHelper : public XRayAtomicGasMix::ScatteringHelper
    {
    public:
        NoScatteringHelper(SimulationItem* /*item*/) {}

        double sectionSca(double /*lambda*/, int /*Z*/) const override { return 0.; }

        void peeloffScattering(double& /*I*/, double& /*lambda*/, int /*Z*/, Direction /*bfk*/,
                               Direction /*bfkobs*/) const override
        {}

        Direction performScattering(double& /*lambda*/, int /*Z*/, Direction /*bfk*/) const override
        {
            return Direction();
        }
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
    class FreeComptonHelper : public XRayAtomicGasMix::ScatteringHelper
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
                _cpf.peeloffScattering(I, lambda, bfk, bfkobs);
            }
            else
            {
                double Q, U, V;
                _dpf.peeloffScattering(I, Q, U, V, bfk, bfkobs, Direction(), nullptr);
            }
        }

        Direction performScattering(double& lambda, int /*Z*/, Direction bfk) const override
        {
            return lambda < comptonWL ? _cpf.performScattering(lambda, bfk) : _dpf.performScattering(bfk, nullptr);
        }
    };
}

////////////////////////////////////////////////////////////////////

// ---- bound-electron Compton scattering helper ----

namespace
{
    // returns the inverse Compton factor for a given scaled energy and scattering angle cosine
    constexpr double inverseComptonfactor(double x, double costheta) { return 1 + x * (1 - costheta); }

    // returns the Compton factor for a given scaled energy and scattering angle cosine
    constexpr double comptonFactor(double x, double costheta) { return 1. / inverseComptonfactor(x, costheta); }

    // returns the value interpolated from the specified table as a function of the momentum transfer parameter
    // q = (E/12.4 keV) sin(theta/2), given the scaled energy x and the sine;
    // logarithmic interpolation is used except for q values near zero
    double interpolateQ(double x, double sintheta2, const Array& qv, const Array& fv)
    {
        double q = scaledEnergyTo12keV * x * sintheta2;
        if (q < 1e-3) return NR::clampedValue<NR::interpolateLinLin>(q, qv, fv);
        return NR::clampedValue<NR::interpolateLogLog>(q, qv, fv);
    }

    // this helper implements bound-electron Compton scattering
    class BoundComptonHelper : public XRayAtomicGasMix::ScatteringHelper
    {
    private:
        vector<Array> _CSv;  // 0: E (keV->1); 1-30: bound Compton cross sections (cm2->m2)
        vector<Array> _SFv;  // 0: q (1); 1-30: incoherent scattering functions (1)
        Random* _random{nullptr};

        // precalculated discretizations
        Array _costhetav = Array(numTheta);
        Array _sinthetav = Array(numTheta);
        Array _sin2thetav = Array(numTheta);
        Array _sintheta2v = Array(numTheta);

    public:
        BoundComptonHelper(SimulationItem* item)
        {
            // load bound Compton cross sections
            _CSv = loadColumns(numAtoms + 1, item, "XRay_CS.txt", "bound Compton data");
            _CSv[0] *= keVtoScaledEnergy;                            // convert from keV to 1
            for (size_t Z = 1; Z <= numAtoms; ++Z) _CSv[Z] *= 1e-4;  // convert from cm2 to m2

            // load incoherent scattering functions
            _SFv = loadColumns(numAtoms + 1, item, "XRay_SF.txt", "bound Compton data");

            // cache random nr generator
            _random = item->find<Random>();

            // construct a theta grid and precalculate values used in generateCosineFromPhaseFunction()
            // to accelerate construction of the cumulative phase function distribution
            for (size_t t = 0; t != numTheta; ++t)
            {
                double theta = t * deltaTheta;
                _costhetav[t] = cos(theta);
                _sinthetav[t] = sin(theta);
                _sin2thetav[t] = _sinthetav[t] * _sinthetav[t];
                _sintheta2v[t] = sin(0.5 * theta);
            }
        }

        double sectionSca(double lambda, int Z) const override
        {
            // interpolate from table, and:
            // - below lower table limit: cross section must be zero so don't clamp values
            // - above upper table limit: does not matter because this limit coincides with the global upper limit
            return NR::value<NR::interpolateLogLog>(scaledEnergy(lambda), _CSv[0], _CSv[Z]);
        }

    private:
        double phaseFunctionValue(double x, double costheta, int Z) const
        {
            constexpr double norm = 3. / 4. * Constants::sigmaThomson();
            double C = comptonFactor(x, costheta);
            double sin2theta = (1 - costheta) * (1 + costheta);
            double phase = C * C * C + C - C * C * sin2theta;
            double section = NR::value<NR::interpolateLogLog>(x, _CSv[0], _CSv[Z]);
            double sintheta2 = sqrt(0.5 * (1 - costheta));
            double incoherent = interpolateQ(x, sintheta2, _SFv[0], _SFv[Z]);
            return norm / section * phase * incoherent;
        }

        double generateCosineFromPhaseFunction(double x, double Z) const
        {
            // construct the normalized cumulative phase function distribution for this x
            Array thetaXv;
            NR::cdf(thetaXv, maxTheta, [this, x, Z](int t) {
                t += 1;
                double C = comptonFactor(x, _costhetav[t]);
                double phase = C * C * C + C - C * C * _sin2thetav[t];
                double incoherent = interpolateQ(x, _sintheta2v[t], _SFv[0], _SFv[Z]);
                return phase * incoherent * _sinthetav[t];
            });

            // draw a random cosine from this distribution
            return _random->cdfLinLin(_costhetav, thetaXv);
        }

    public:
        void peeloffScattering(double& I, double& lambda, int Z, Direction bfk, Direction bfkobs) const override
        {
            double x = scaledEnergy(lambda);

            // calculate the value of the phase function
            double costheta = Vec::dot(bfk, bfkobs);
            double value = phaseFunctionValue(x, costheta, Z);

            // accumulate the weighted sum in the intensity
            I += value;

            // adjust the wavelength
            lambda *= inverseComptonfactor(x, costheta);
        }

        Direction performScattering(double& lambda, int Z, Direction bfk) const override
        {
            double x = scaledEnergy(lambda);

            // sample a scattering angle from the phase function
            double costheta = generateCosineFromPhaseFunction(x, Z);

            // adjust the wavelength
            lambda *= inverseComptonfactor(x, costheta);

            // determine the new propagation direction
            return _random->direction(bfk, costheta);
        }
    };
}

////////////////////////////////////////////////////////////////////

// ---- smooth Rayleigh scattering helper ----

namespace
{
    // this helper implements smooth Rayleigh scattering;
    // below the energy limit of the tabulated data, use Thomson scattering instead
    class SmoothRayleighHelper : public XRayAtomicGasMix::ScatteringHelper
    {
    private:
        vector<Array> _RSSv;  // 0: E (keV->1); 1-30: smooth Rayleigh cross sections (cm2->m2)
        vector<Array> _FFv;   // 0: q (1); 1-30: atomic form factors (1)
        Random* _random{nullptr};
        DipolePhaseFunction _dpf;

        // precalculated discretizations
        Array _costhetav = Array(numTheta);
        Array _cos2thetav = Array(numTheta);
        Array _sinthetav = Array(numTheta);
        Array _sintheta2v = Array(numTheta);

    public:
        SmoothRayleighHelper(SimulationItem* item)
        {
            // load smooth Rayleigh cross sections
            _RSSv = loadColumns(numAtoms + 1, item, "XRay_RSS.txt", "smooth Rayleigh data");
            _RSSv[0] *= keVtoScaledEnergy;                            // convert from keV to 1
            for (size_t Z = 1; Z <= numAtoms; ++Z) _RSSv[Z] *= 1e-4;  // convert from cm2 to m2

            // load atomic form factors
            _FFv = loadColumns(numAtoms + 1, item, "XRay_FF.txt", "smooth Rayleigh data");

            // cache random nr generator and initialize the Thomson helper
            _random = item->find<Random>();
            _dpf.initialize(_random);

            // construct a theta grid and precalculate values used in generateCosineFromPhaseFunction()
            // to accelerate construction of the cumulative phase function distribution
            for (size_t t = 0; t != numTheta; ++t)
            {
                double theta = t * deltaTheta;
                _costhetav[t] = cos(theta);
                _cos2thetav[t] = _costhetav[t] * _costhetav[t];
                _sinthetav[t] = sin(theta);
                _sintheta2v[t] = sin(0.5 * theta);
            }
        }

        double sectionSca(double lambda, int Z) const override
        {
            // interpolate from table, and:
            // - below lower table limit: use Z^2 * Thomson scattering
            // - above upper table limit: does not matter because this limit coincides with the global upper limit
            double x = scaledEnergy(lambda);
            if (x < _RSSv[0][0]) return Z * Z * Constants::sigmaThomson();
            return NR::value<NR::interpolateLogLog>(x, _RSSv[0], _RSSv[Z]);
        }

    private:
        double phaseFunctionValue(double x, double costheta, int Z) const
        {
            constexpr double norm = 3. / 4. * Constants::sigmaThomson();
            double phase = 1. + costheta * costheta;
            double section = NR::value<NR::interpolateLogLog>(x, _RSSv[0], _RSSv[Z]);
            double sintheta2 = sqrt(0.5 * (1 - costheta));
            double form = interpolateQ(x, sintheta2, _FFv[0], _FFv[Z]);
            return norm / section * phase * form * form;
        }

        double generateCosineFromPhaseFunction(double x, double Z) const
        {
            // construct the normalized cumulative phase function distribution for this x
            Array thetaXv;
            NR::cdf(thetaXv, maxTheta, [this, x, Z](int t) {
                t += 1;
                double phase = 1. + _cos2thetav[t];
                double form = interpolateQ(x, _sintheta2v[t], _FFv[0], _FFv[Z]);
                return phase * form * form * _sinthetav[t];
            });

            // draw a random cosine from this distribution
            return _random->cdfLinLin(_costhetav, thetaXv);
        }

    public:
        void peeloffScattering(double& I, double& lambda, int Z, Direction bfk, Direction bfkobs) const override
        {
            double x = scaledEnergy(lambda);

            // for low energies use Thomson scattering
            if (x < _RSSv[0][0])
            {
                double Q, U, V;
                _dpf.peeloffScattering(I, Q, U, V, bfk, bfkobs, Direction(), nullptr);
            }

            // otherwise use Rayleigh scattering
            else
            {
                // calculate the value of the phase function
                double costheta = Vec::dot(bfk, bfkobs);
                double value = phaseFunctionValue(x, costheta, Z);

                // accumulate the weighted sum in the intensity
                I += value;
            }
        }

        Direction performScattering(double& lambda, int Z, Direction bfk) const override
        {
            double x = scaledEnergy(lambda);

            // for low energies use Thomson scattering
            if (x < _RSSv[0][0]) return _dpf.performScattering(bfk, nullptr);

            // otherwise use Rayleigh scattering
            return _random->direction(bfk, generateCosineFromPhaseFunction(x, Z));
        }
    };
}

////////////////////////////////////////////////////////////////////

// ---- anomalous Rayleigh scattering helper ----

namespace
{
    // this helper implements anomalous Rayleigh scattering
    // below the energy limit of the tabulated data, use Thomson scattering instead
    class AnomalousRayleighHelper : public XRayAtomicGasMix::ScatteringHelper
    {
    private:
        vector<Array> _RSAv;  // 2*Z: E (keV->1); 2*Z+1: anomalous Rayleigh cross sections (cm2->m2)
        vector<Array> _FFv;   // 0: q (1); 1-30: atomic form factors (1)
        vector<Array> _F1v;   // 2*Z: E (keV->1); 2*Z+1: Real anomalous scattering function (1)
        vector<Array> _F2v;   // 2*Z: E (keV->1); 2*Z+1: Imaginary anomalous scattering function (1)
        Random* _random{nullptr};
        DipolePhaseFunction _dpf;

        // precalculated discretizations
        Array _costhetav = Array(numTheta);
        Array _cos2thetav = Array(numTheta);
        Array _sinthetav = Array(numTheta);
        Array _sintheta2v = Array(numTheta);

    public:
        AnomalousRayleighHelper(SimulationItem* item)
        {
            // load anomalous Rayleigh cross sections, atomic form factors and anomalous scattering functions
            _RSAv = loadColumns(2 * numAtoms + 2, item, "XRay_RSA.txt", "anomalous Rayleigh data");
            _FFv = loadColumns(numAtoms + 1, item, "XRay_FF.txt", "anomalous Rayleigh data");
            _F1v = loadColumns(2 * numAtoms + 2, item, "XRay_F1.txt", "anomalous Rayleigh data");
            _F2v = loadColumns(2 * numAtoms + 2, item, "XRay_F2.txt", "anomalous Rayleigh data");

            // convert units
            for (size_t Z = 1; Z <= numAtoms; ++Z)
            {
                _RSAv[2 * Z] *= keVtoScaledEnergy;  // convert from keV to 1
                _F1v[2 * Z] *= keVtoScaledEnergy;   // convert from keV to 1
                _F2v[2 * Z] *= keVtoScaledEnergy;   // convert from keV to 1
                _RSAv[2 * Z + 1] *= 1e-4;           // convert from cm2 to m2
            }

            // cache random nr generator and initialize the Thomson helper
            _random = item->find<Random>();
            _dpf.initialize(_random);

            // construct a theta grid and precalculate values used in generateCosineFromPhaseFunction()
            // to accelerate construction of the cumulative phase function distribution
            for (size_t t = 0; t != numTheta; ++t)
            {
                double theta = t * deltaTheta;
                _costhetav[t] = cos(theta);
                _cos2thetav[t] = _costhetav[t] * _costhetav[t];
                _sinthetav[t] = sin(theta);
                _sintheta2v[t] = sin(0.5 * theta);
            }
        }

        double sectionSca(double lambda, int Z) const override
        {
            // interpolate from table, and:
            // - below lower table limit: use Z^2 * Thomson scattering
            // - above upper table limit: use clamped value
            double x = scaledEnergy(lambda);
            if (x < _RSAv[2 * Z][0]) return Z * Z * Constants::sigmaThomson();
            return NR::clampedValue<NR::interpolateLogLog>(x, _RSAv[2 * Z], _RSAv[2 * Z + 1]);
        }

    private:
        double phaseFunctionValue(double x, double costheta, int Z) const
        {
            constexpr double norm = 3. / 4. * Constants::sigmaThomson();
            double phase = 1. + costheta * costheta;
            double section = NR::clampedValue<NR::interpolateLogLog>(x, _RSAv[2 * Z], _RSAv[2 * Z + 1]);
            double sintheta2 = sqrt(0.5 * (1 - costheta));
            double form = interpolateQ(x, sintheta2, _FFv[0], _FFv[Z]);
            double form1 = NR::clampedValue<NR::interpolateLogLin>(x, _F1v[2 * Z], _F1v[2 * Z + 1]);  // negative values
            double form2 = NR::clampedValue<NR::interpolateLogLog>(x, _F2v[2 * Z], _F2v[2 * Z + 1]);
            double formsum = form + form1;
            return norm / section * phase * (formsum * formsum + form2 * form2);
        }

        double generateCosineFromPhaseFunction(double x, double Z) const
        {
            // construct the normalized cumulative phase function distribution for this x
            Array thetaXv;
            NR::cdf(thetaXv, maxTheta, [this, x, Z](int t) {
                t += 1;
                double phase = 1. + _cos2thetav[t];
                double form = interpolateQ(x, _sintheta2v[t], _FFv[0], _FFv[Z]);
                double form1 = NR::clampedValue<NR::interpolateLogLin>(x, _F1v[2 * Z], _F1v[2 * Z + 1]);
                double form2 = NR::clampedValue<NR::interpolateLogLog>(x, _F2v[2 * Z], _F2v[2 * Z + 1]);
                double formsum = form + form1;
                return phase * (formsum * formsum + form2 * form2) * _sinthetav[t];
            });

            // draw a random cosine from this distribution
            return _random->cdfLinLin(_costhetav, thetaXv);
        }

    public:
        void peeloffScattering(double& I, double& lambda, int Z, Direction bfk, Direction bfkobs) const override
        {
            double x = scaledEnergy(lambda);

            // for low energies use Thomson scattering
            if (x < _RSAv[2 * Z][0])
            {
                double Q, U, V;
                _dpf.peeloffScattering(I, Q, U, V, bfk, bfkobs, Direction(), nullptr);
            }

            // otherwise use Rayleigh scattering
            else
            {
                // calculate the value of the phase function
                double costheta = Vec::dot(bfk, bfkobs);
                double value = phaseFunctionValue(x, costheta, Z);

                // accumulate the weighted sum in the intensity
                I += value;
            }
        }

        Direction performScattering(double& lambda, int Z, Direction bfk) const override
        {
            double x = scaledEnergy(lambda);

            // for low energies use Thomson scattering
            if (x < _RSAv[2 * Z][0]) return _dpf.performScattering(bfk, nullptr);

            // otherwise use Rayleigh scattering
            return _random->direction(bfk, generateCosineFromPhaseFunction(x, Z));
        }
    };
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
    auto fluorescenceParams = loadStruct<FluorescenceParams, 4>(this, "XRay_FL.txt", "fluorescence data");

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
            sigma += (_ray->sectionSca(lambda, Z) + _com->sectionSca(lambda, Z)) * atomv[Z - 1].abund;
        }

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

    // calculate and store the fluorescence emission wavelengths
    _lambdafluov.reserve(fluorescenceParams.size());
    for (const auto& params : fluorescenceParams) _lambdafluov.push_back(wavelengthToFromEnergy(params.E));

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
            contribv[Z - 1] = _ray->sectionSca(lambda, Z) * atomv[Z - 1].abund;
            contribv[numAtoms + Z - 1] = _com->sectionSca(lambda, Z) * atomv[Z - 1].abund;
        }

        // fluorescence: iterate over both cross section and fluorescence parameter sets in sync
        auto flp = fluorescenceParams.begin();
        int index = 0;
        for (const auto& csp : crossSectionParams)
        {
            auto sigmoid = sigmoidv[index++];

            // process all fluorescence parameter sets matching this cross section set
            while (flp != fluorescenceParams.end() && flp->Z == csp.Z && flp->n == csp.n)
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

bool XRayAtomicGasMix::hasScatteringDispersion() const
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
    }
}

////////////////////////////////////////////////////////////////////

void XRayAtomicGasMix::peeloffScattering(double& I, double& /*Q*/, double& /*U*/, double& /*V*/, double& lambda,
                                         Direction bfkobs, Direction /*bfky*/, const MaterialState* /*state*/,
                                         const PhotonPacket* pp) const
{
    // draw a random scattering channel and atom velocity, unless a previous peel-off stored this already
    auto scatinfo = const_cast<PhotonPacket*>(pp)->getScatteringInfo();
    setScatteringInfoIfNeeded(scatinfo, lambda);

    // if we have dispersion, for electron scattering, adjust the incoming wavelength to the electron rest frame
    if (temperature() > 0. && scatinfo->species < static_cast<int>(2 * numAtoms))
        lambda = PhotonPacket::shiftedReceptionWavelength(lambda, pp->direction(), scatinfo->velocity);

    // Rayleigh scattering in electron rest frame
    if (scatinfo->species < static_cast<int>(numAtoms))
    {
        _ray->peeloffScattering(I, lambda, scatinfo->species + 1, pp->direction(), bfkobs);
    }

    // Compton scattering in electron rest frame
    else if (scatinfo->species < static_cast<int>(2 * numAtoms))
    {
        _com->peeloffScattering(I, lambda, scatinfo->species - numAtoms + 1, pp->direction(), bfkobs);
    }

    // fluorescence
    else
    {
        // isotropic emission, so the bias weight is trivially 1
        I += 1.;

        // update the photon packet wavelength to the wavelength of this fluorescence transition
        lambda = _lambdafluov[scatinfo->species - 2 * numAtoms];
    }

    // if we have dispersion, Doppler-shift the outgoing wavelength from the electron rest frame
    if (temperature() > 0.) lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfkobs, scatinfo->velocity);
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

    // Rayleigh scattering: determine the new propagation direction
    if (scatinfo->species < static_cast<int>(numAtoms))
    {
        bfknew = _ray->performScattering(lambda, scatinfo->species + 1, pp->direction());
    }

    // Compton scattering: determine the new propagation direction and wavelength
    else if (scatinfo->species < static_cast<int>(2 * numAtoms))
    {
        bfknew = _com->performScattering(lambda, scatinfo->species - numAtoms + 1, pp->direction());
    }

    // fluorescence
    else
    {
        // update the photon packet wavelength to the wavelength of this fluorescence transition
        lambda = _lambdafluov[scatinfo->species - 2 * numAtoms];

        // draw a random, isotropic outgoing direction
        bfknew = random()->direction();
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
