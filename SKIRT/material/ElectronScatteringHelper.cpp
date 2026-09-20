/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ElectronScatteringHelper.hpp"
#include "Constants.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "TextInFile.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // ---- helper functions ----

    // convert photon energy in eV to and from wavelength in m (same conversion in both directions)
    constexpr double wavelengthToFromEnergy(double x)
    {
        constexpr double front = Constants::h() * Constants::c() / Constants::Qelectron();
        return front / x;
    }

    // convert wavelength to scaled photon energy: h nu / (m_e c^2)
    constexpr double scaledEnergy(double lambda)
    {
        constexpr double front = Constants::h() / Constants::Melectron() / Constants::c();
        return front / lambda;
    }

    // returns the inverse Compton factor for a given scaled energy and scattering angle cosine
    constexpr double inverseComptonFactor(double x, double costheta)
    {
        return 1. + x * (1. - costheta);
    }

    // returns the Compton factor for a given scaled energy and scattering angle cosine
    constexpr double comptonFactor(double x, double costheta)
    {
        return 1. / inverseComptonFactor(x, costheta);
    }

    // returns the logarithms of the ratios of successive values in each of the columns of the specified table,
    // as a vector with the same layout as the table but with one value less in each array;
    // these are the terms of NR::interpolateLogLog() that depend on the table only and not on the point being
    // interpolated, so it saves time to calculate them once when the table is loaded
    vector<Array> logRatios(const vector<Array>& columns)
    {
        vector<Array> result;
        result.reserve(columns.size());
        for (const Array& column : columns)
        {
            Array ratios(column.size() ? column.size() - 1 : 0);
            for (size_t i = 0; i != ratios.size(); ++i) ratios[i] = log(column[i + 1] / column[i]);
            result.push_back(ratios);
        }
        return result;
    }

    // multiplicator to convert scaled energy to energy in units of 12.4 keV
    constexpr double scaledEnergyTo12keV =
        (Constants::Melectron() * Constants::c() * Constants::c()) / (12.4e3 * Constants::Qelectron());

    // returns the value at q interpolated between the points i and i+1 of a table with abscissae qv and function values fv,
    // given the precalculated logarithms of the ratios of successive values in qv and fv (see logRatios());
    // logarithmic interpolation is used except for q values near zero, with the same result as NR::interpolateLogLog()
    double interpolateInterval(double q, size_t i, const Array& qv, const Array& fv, const Array& lqv, const Array& lfv)
    {
        if (q < 1e-3) return NR::interpolateLinLin(q, qv[i], qv[i + 1], fv[i], fv[i + 1]);
        double f1 = fv[i];
        double f2 = fv[i + 1];
        if (f1 <= 0 || f2 <= 0)
        {
            if (q == qv[i]) return f1;
            if (q == qv[i + 1]) return f2;
            return 0.;
        }
        return f1 * exp(log(q / qv[i]) / lqv[i] * lfv[i]);
    }

    // returns the value interpolated from the specified table as a function of the momentum transfer parameter
    // q = (E/12.4 keV) sin(theta/2), given the scaled energy x and the sine;
    // the table contains q in column 0 and the function values for atomic number Z in column Z, and logRatiov contains
    // the precalculated logarithms of the ratios of successive values in each of these columns (see logRatios())
    double interpolateQ(double x, double sintheta2, const vector<Array>& table, const vector<Array>& logRatiov, int Z)
    {
        const Array& qv = table[0];
        const Array& fv = table[Z];
        double q = scaledEnergyTo12keV * x * sintheta2;

        // clamp outside of the table and locate the interval inside, as NR::clampedValue() does
        size_t n = qv.size();
        if (q < qv[0]) return fv[0];
        if (qv[n - 1] < q) return fv[n - 1];
        size_t i = static_cast<size_t>(NR::locate(qv, q));
        return interpolateInterval(q, i, qv, fv, logRatiov[0], logRatiov[Z]);
    }

    // returns the values interpolated from the specified table for all sines of half the scattering angle in sintheta2v,
    // given the scaled energy x; see above;
    // this is much faster than interpolating the values one by one if the sines are ascending, as they are in the loops
    // over the scattering angle, because the interval containing q is then found by walking forward through the table
    // from the interval found for the previous value, rather than by a binary search for each value (the results do not
    // depend on the order of the values)
    Array interpolateQ(double x, const Array& sintheta2v, const vector<Array>& table, const vector<Array>& logRatiov,
                       int Z)
    {
        const Array& qv = table[0];
        const Array& fv = table[Z];
        size_t n = qv.size();

        Array resultv(sintheta2v.size());
        size_t i = 0;
        bool found = false;
        for (size_t k = 0; k != sintheta2v.size(); ++k)
        {
            double q = scaledEnergyTo12keV * x * sintheta2v[k];

            // clamp outside of the table, as above
            if (q < qv[0])
                resultv[k] = fv[0];
            else if (qv[n - 1] < q)
                resultv[k] = fv[n - 1];
            else
            {
                // search for the interval containing q, or walk forward from the interval found for the previous value
                if (!found || q < qv[i])
                {
                    i = static_cast<size_t>(NR::locate(qv, q));
                    found = true;
                }
                else
                    while (i + 2 < n && qv[i + 1] <= q) ++i;
                resultv[k] = interpolateInterval(q, i, qv, fv, logRatiov[0], logRatiov[Z]);
            }
        }
        return resultv;
    }

    // number of supported atoms; the data provided in the resource files must match this number
    constexpr size_t numAtoms = 30;

    // transition wavelength from Compton to Thomson scattering
    constexpr double comptonWL = wavelengthToFromEnergy(100.);  // 0.1 keV or 12.4 nm

    // multiplicator to convert energy in keV to scaled energy E / (m_e c^2)
    constexpr double keVtoScaledEnergy =
        (1e3 * Constants::Qelectron()) / (Constants::Melectron() * Constants::c() * Constants::c());

    // discretization of the phase function over scattering angle: theta from 0 to pi, index t
    constexpr size_t numTheta = 361;
    constexpr size_t maxTheta = numTheta - 1;
    constexpr double deltaTheta = M_PI / maxTheta;

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

namespace ElectronScatteringHelper
{
    // ---- base class for scattering helpers ----

    Helper::~Helper() {}

    // peel-off unpolarized scattering event: this in helpers that don't support polarization
    void Helper::peeloffScattering(double& /*I*/, double& /*lambda*/, int /*Z*/, int /*N*/, Direction /*bfk*/,
                                   Direction /*bfkobs*/) const
    {
        // default implementation does nothing
    }

    ////////////////////////////////////////////////////////////////////

    // perform unpolarized scattering event: this in helpers that don't support polarization
    Direction Helper::performScattering(double& /*lambda*/, int /*Z*/, int /*N*/, Direction /*bfk*/) const
    {
        // default implementation returns null vector
        return Direction();
    }

    ////////////////////////////////////////////////////////////////////

    // peel-off polarized scattering event: this in helpers that do support polarization
    void Helper::peeloffScattering(double& I, double& /*Q*/, double& /*U*/, double& /*V*/, double& lambda, int Z, int N,
                                   Direction bfk, Direction bfkobs, Direction /*bfky*/,
                                   const StokesVector* /*sv*/) const
    {
        // default implementation calls unpolarized version
        peeloffScattering(I, lambda, Z, N, bfk, bfkobs);
    }

    ////////////////////////////////////////////////////////////////////

    // perform polarized scattering event: this in helpers that do support polarization
    Direction Helper::performScattering(double& lambda, int Z, int N, Direction bfk, StokesVector* /*sv*/) const
    {
        // default implementation calls unpolarized version
        return performScattering(lambda, Z, N, bfk);
    }

    ////////////////////////////////////////////////////////////////////

    // ---- no scattering helper ----

    NoScatteringHelper::NoScatteringHelper(SimulationItem* /*item*/) {}

    ////////////////////////////////////////////////////////////////////

    double NoScatteringHelper::sectionSca(double /*lambda*/, int /*Z*/, int /*N*/) const
    {
        return 0.;
    }

    ////////////////////////////////////////////////////////////////////

    // ---- free-electron Compton scattering helper ----

    FreeComptonHelper::FreeComptonHelper(SimulationItem* item)
    {
        auto random = item->find<Random>();
        _cpf.initialize(random);
        _dpf.initialize(random);
    }

    ////////////////////////////////////////////////////////////////////

    double FreeComptonHelper::sectionSca(double lambda, int Z, int /*N*/) const
    {
        double sigma = Z * Constants::sigmaThomson();
        if (lambda < comptonWL) sigma *= _cpf.sectionSca(lambda);
        return sigma;
    }

    ////////////////////////////////////////////////////////////////////

    void FreeComptonHelper::peeloffScattering(double& I, double& lambda, int /*Z*/, int /*N*/, Direction bfk,
                                              Direction bfkobs) const
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

    ////////////////////////////////////////////////////////////////////

    Direction FreeComptonHelper::performScattering(double& lambda, int /*Z*/, int /*N*/, Direction bfk) const
    {
        return lambda < comptonWL ? _cpf.performScattering(lambda, bfk, nullptr) : _dpf.performScattering(bfk, nullptr);
    }

    ////////////////////////////////////////////////////////////////////

    // ---- free-electron Compton with polarization scattering helper ----

    FreeComptonWithPolarizationHelper::FreeComptonWithPolarizationHelper(SimulationItem* item)
    {
        auto random = item->find<Random>();
        _cpf.initialize(random, true);
        _dpf.initialize(random, true);
    }

    ////////////////////////////////////////////////////////////////////

    double FreeComptonWithPolarizationHelper::sectionSca(double lambda, int Z, int /*N*/) const
    {
        double sigma = Z * Constants::sigmaThomson();
        if (lambda < comptonWL) sigma *= _cpf.sectionSca(lambda);
        return sigma;
    }

    ////////////////////////////////////////////////////////////////////

    void FreeComptonWithPolarizationHelper::peeloffScattering(double& I, double& Q, double& U, double& V,
                                                              double& lambda, int /*Z*/, int /*N*/, Direction bfk,
                                                              Direction bfkobs, Direction bfky,
                                                              const StokesVector* sv) const
    {
        lambda < comptonWL ? _cpf.peeloffScattering(I, Q, U, V, lambda, bfk, bfkobs, bfky, sv)
                           : _dpf.peeloffScattering(I, Q, U, V, bfk, bfkobs, bfky, sv);
    }

    ////////////////////////////////////////////////////////////////////

    Direction FreeComptonWithPolarizationHelper::performScattering(double& lambda, int /*Z*/, int /*N*/, Direction bfk,
                                                                   StokesVector* sv) const
    {
        return lambda < comptonWL ? _cpf.performScattering(lambda, bfk, sv) : _dpf.performScattering(bfk, sv);
    }

    ////////////////////////////////////////////////////////////////////

    // ---- bound-electron Compton scattering helper ----

    BoundComptonHelper::BoundComptonHelper(SimulationItem* item)
        : _costhetav(numTheta), _sinthetav(numTheta), _sin2thetav(numTheta), _sintheta2v(numTheta)
    {
        // load bound Compton cross sections
        _CSv = loadColumns(numAtoms + 1, item, "XRay_CS.txt", "bound Compton data");
        _CSv[0] *= keVtoScaledEnergy;                            // convert from keV to 1
        for (size_t Z = 1; Z <= numAtoms; ++Z) _CSv[Z] *= 1e-4;  // convert from cm2 to m2

        // load incoherent scattering functions and precalculate the logarithms needed to interpolate them
        _SFv = loadColumns(numAtoms + 1, item, "XRay_SF.txt", "bound Compton data");
        _logRatioSFv = logRatios(_SFv);

        // load pdfs for projected momentum of target electron
        _CPv = loadColumns(numAtoms + 1, item, "XRay_CP.txt", "bound Compton data");
        _CPv[0] *= keVtoScaledEnergy;  // convert from keV to 1

        // load ionization energies
        _IBv = loadColumns(1, item, "XRay_IB.txt", "bound Compton data");
        _IBv[0] *= keVtoScaledEnergy;  // convert from keV to 1

        // precalculate cumulative distributions for target electron momentum
        _cumRange.set(_CPv[0][0], _CPv[0][_CPv[0].size() - 1]);
        Array xv, pv, Pv;
        NR::cdf<NR::interpolateLinLin>(xv, pv, Pv, _CPv[0], _CPv[1], _cumRange);
        _cumCPv.push_back(xv);
        _cumCPv.push_back(Pv);
        for (size_t Z = 2; Z <= numAtoms; ++Z)
        {
            NR::cdf<NR::interpolateLinLin>(xv, pv, Pv, _CPv[0], _CPv[Z], _cumRange);
            _cumCPv.push_back(Pv);
        }

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

        // cache random nr generator
        _random = item->find<Random>();
    }

    ////////////////////////////////////////////////////////////////////

    double BoundComptonHelper::sectionSca(double lambda, int Z, int /*N*/) const
    {
        // interpolate from table, and:
        // - below lower table limit: cross section must be zero so don't clamp values
        // - above upper table limit: does not matter because this limit coincides with the global upper limit
        return NR::value<NR::interpolateLogLog>(scaledEnergy(lambda), _CSv[0], _CSv[Z]);
    }

    ////////////////////////////////////////////////////////////////////

    double BoundComptonHelper::phaseFunctionValue(double x, double costheta, int Z) const
    {
        constexpr double norm = 3. / 4. * Constants::sigmaThomson();
        double C = comptonFactor(x, costheta);
        double sin2theta = (1 - costheta) * (1 + costheta);
        double phase = C * C * C + C - C * C * sin2theta;
        double section = NR::value<NR::interpolateLogLog>(x, _CSv[0], _CSv[Z]);
        double sintheta2 = sqrt(0.5 * (1 - costheta));
        double incoherent = interpolateQ(x, sintheta2, _SFv, _logRatioSFv, Z);
        return norm / section * phase * incoherent;
    }

    ////////////////////////////////////////////////////////////////////

    double BoundComptonHelper::generateCosineFromPhaseFunction(double x, double Z) const
    {
        // construct the normalized cumulative phase function distribution for this x
        Array thetaXv;
        Array incoherentv = interpolateQ(x, _sintheta2v, _SFv, _logRatioSFv, Z);
        NR::cdf(thetaXv, maxTheta, [this, x, &incoherentv](int t) {
            t += 1;
            double C = comptonFactor(x, _costhetav[t]);
            double phase = C * C * C + C - C * C * _sin2thetav[t];
            double incoherent = incoherentv[t];
            return phase * incoherent * _sinthetav[t];
        });

        // draw a random cosine from this distribution
        return _random->cdfLinLin(_costhetav, thetaXv);
    }

    ////////////////////////////////////////////////////////////////////

    // sample a target electron momentum from the distribution with the given maximum
    double BoundComptonHelper::sampleMomentum(double pmax, double Z) const
    {
        // maximum momentum is below the range of the tabulated pdf -> simply return the maximum momentum
        // (we estimate that this happens for less than 0.1 % of the events)
        if (pmax <= _cumRange.min()) return pmax;

        // maximum momentum is on the left side of the peak in the tabulated pdf;
        // using the rejection technique on the full-range pdf is very inefficient
        // because the majority of the generated samples would be rejected
        // --> reconstruct a cumulative pdf with the appropriate range and use numerical inversion
        // (we estimate that this happens for less than 10% of the events)
        if (pmax <= _cumRange.mid())
        {
            Array xv, pv, Pv;
            NR::cdf<NR::interpolateLinLin>(xv, pv, Pv, _CPv[0], _CPv[Z], Range(_cumRange.min(), pmax));
            return _random->cdfLinLin(xv, Pv);
        }

        // maximum momentum is on the right side of the peak in the tabulated pdf, possibly even out of range;
        // using the rejection technique on top of numerical inversion for the full-range pdf now is efficient
        // and quite fast because we can use the precalculated cumulative pdf
        while (true)
        {
            double p = _random->cdfLinLin(_cumCPv[0], _cumCPv[Z]);
            if (p <= pmax) return p;
        }
    }

    ////////////////////////////////////////////////////////////////////

    // returns the augmented inverse Compton factor
    double BoundComptonHelper::augmentedInverseComptonFactor(double x, double costheta, double Z) const
    {
        // precalculate some values
        double costheta1 = 1. - costheta;
        double sintheta22 = 2. * sqrt(0.5 * costheta1);  // twice the half-angle sine

        // calculate the maximum target electron momentum (in scaled energy units)
        double b = _IBv[0][Z - 1];  // scaled ionization energy
        double xminb = (x - b);
        double pmax = (x * xminb * costheta1 - b) / (xminb * sintheta22);

        // sample a target electron momentum from the distribution with the given maximum
        double p = sampleMomentum(pmax, Z);

        // calculate the augmented inverse Compton factor
        return 1. + x * costheta1 - p * sintheta22;
    }

    ////////////////////////////////////////////////////////////////////

    void BoundComptonHelper::peeloffScattering(double& I, double& lambda, int Z, int /*N*/, Direction bfk,
                                               Direction bfkobs) const
    {
        double x = scaledEnergy(lambda);

        // calculate the value of the phase function
        double costheta = Vec::dot(bfk, bfkobs);
        double value = phaseFunctionValue(x, costheta, Z);

        // accumulate the weighted sum in the intensity
        I += value;

        // adjust the wavelength
        lambda *= augmentedInverseComptonFactor(x, costheta, Z);
    }

    ////////////////////////////////////////////////////////////////////

    Direction BoundComptonHelper::performScattering(double& lambda, int Z, int /*N*/, Direction bfk) const
    {
        double x = scaledEnergy(lambda);

        // sample a scattering angle from the phase function
        double costheta = generateCosineFromPhaseFunction(x, Z);

        // adjust the wavelength
        lambda *= augmentedInverseComptonFactor(x, costheta, Z);

        // determine the new propagation direction
        return _random->direction(bfk, costheta);
    }

    ////////////////////////////////////////////////////////////////////

    // ---- free-bound Compton scattering helper ----

    FreeBoundComptonHelper::FreeBoundComptonHelper(SimulationItem* item) : _free(item), _bound(item)
    {
        _random = item->find<Random>();
    }

    ////////////////////////////////////////////////////////////////////

    double FreeBoundComptonHelper::sectionSca(double lambda, int Z, int N) const
    {
        double b = static_cast<double>(N) / Z;
        return (1. - b) * _free.sectionSca(lambda, Z, N) + b * _bound.sectionSca(lambda, Z, N);
    }

    ////////////////////////////////////////////////////////////////////

    void FreeBoundComptonHelper::peeloffScattering(double& I, double& lambda, int Z, int N, Direction bfk,
                                                   Direction bfkobs) const
    {
        double b = static_cast<double>(N) / Z;
        double sigmaFree = (1. - b) * _free.sectionSca(lambda, Z, N);
        double sigmaBound = b * _bound.sectionSca(lambda, Z, N);

        double pFree = sigmaFree / (sigmaFree + sigmaBound);
        if (pFree < _random->uniform())
            return _bound.peeloffScattering(I, lambda, Z, N, bfk, bfkobs);
        else
            return _free.peeloffScattering(I, lambda, Z, N, bfk, bfkobs);
    }

    ////////////////////////////////////////////////////////////////////

    Direction FreeBoundComptonHelper::performScattering(double& lambda, int Z, int N, Direction bfk) const
    {
        double b = static_cast<double>(N) / Z;
        double sigmaFree = (1. - b) * _free.sectionSca(lambda, Z, N);
        double sigmaBound = b * _bound.sectionSca(lambda, Z, N);

        double pFree = sigmaFree / (sigmaFree + sigmaBound);
        if (pFree < _random->uniform())
            return _bound.performScattering(lambda, Z, N, bfk);
        else
            return _free.performScattering(lambda, Z, N, bfk);
    }

    ////////////////////////////////////////////////////////////////////

    // ---- smooth Rayleigh scattering helper ----

    SmoothRayleighHelper::SmoothRayleighHelper(SimulationItem* item)
        : _costhetav(numTheta), _cos2thetav(numTheta), _sinthetav(numTheta), _sintheta2v(numTheta)
    {
        // load smooth Rayleigh cross sections
        _RSSv = loadColumns(numAtoms + 1, item, "XRay_RSS.txt", "smooth Rayleigh data");
        _RSSv[0] *= keVtoScaledEnergy;                            // convert from keV to 1
        for (size_t Z = 1; Z <= numAtoms; ++Z) _RSSv[Z] *= 1e-4;  // convert from cm2 to m2

        // load atomic form factors and precalculate the logarithms needed to interpolate them
        _FFv = loadColumns(numAtoms + 1, item, "XRay_FF.txt", "smooth Rayleigh data");
        _logRatioFFv = logRatios(_FFv);

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

    ////////////////////////////////////////////////////////////////////

    double SmoothRayleighHelper::sectionSca(double lambda, int Z, int /*N*/) const
    {
        // interpolate from table, and:
        // - below lower table limit: use Z^2 * Thomson scattering
        // - above upper table limit: does not matter because this limit coincides with the global upper limit
        double x = scaledEnergy(lambda);
        if (x < _RSSv[0][0]) return Z * Z * Constants::sigmaThomson();
        return NR::value<NR::interpolateLogLog>(x, _RSSv[0], _RSSv[Z]);
    }

    ////////////////////////////////////////////////////////////////////

    double SmoothRayleighHelper::phaseFunctionValue(double x, double costheta, int Z) const
    {
        constexpr double norm = 3. / 4. * Constants::sigmaThomson();
        double phase = 1. + costheta * costheta;
        double section = NR::value<NR::interpolateLogLog>(x, _RSSv[0], _RSSv[Z]);
        double sintheta2 = sqrt(0.5 * (1 - costheta));
        double form = interpolateQ(x, sintheta2, _FFv, _logRatioFFv, Z);
        return norm / section * phase * form * form;
    }

    ////////////////////////////////////////////////////////////////////

    double SmoothRayleighHelper::generateCosineFromPhaseFunction(double x, double Z) const
    {
        // construct the normalized cumulative phase function distribution for this x
        Array thetaXv;
        Array formv = interpolateQ(x, _sintheta2v, _FFv, _logRatioFFv, Z);
        NR::cdf(thetaXv, maxTheta, [this, &formv](int t) {
            t += 1;
            double phase = 1. + _cos2thetav[t];
            double form = formv[t];
            return phase * form * form * _sinthetav[t];
        });

        // draw a random cosine from this distribution
        return _random->cdfLinLin(_costhetav, thetaXv);
    }

    ////////////////////////////////////////////////////////////////////

    void SmoothRayleighHelper::peeloffScattering(double& I, double& lambda, int Z, int /*N*/, Direction bfk,
                                                 Direction bfkobs) const
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

    ////////////////////////////////////////////////////////////////////

    Direction SmoothRayleighHelper::performScattering(double& lambda, int Z, int /*N*/, Direction bfk) const
    {
        double x = scaledEnergy(lambda);

        // for low energies use Thomson scattering
        if (x < _RSSv[0][0]) return _dpf.performScattering(bfk, nullptr);

        // otherwise use Rayleigh scattering
        return _random->direction(bfk, generateCosineFromPhaseFunction(x, Z));
    }

    ////////////////////////////////////////////////////////////////////

    // ---- anomalous Rayleigh scattering helper ----

    AnomalousRayleighHelper::AnomalousRayleighHelper(SimulationItem* item)
        : _costhetav(numTheta), _cos2thetav(numTheta), _sinthetav(numTheta), _sintheta2v(numTheta)
    {
        // load anomalous Rayleigh cross sections, atomic form factors and anomalous scattering functions
        _RSAv = loadColumns(2 * numAtoms + 2, item, "XRay_RSA.txt", "anomalous Rayleigh data");
        _FFv = loadColumns(numAtoms + 1, item, "XRay_FF.txt", "anomalous Rayleigh data");
        _F1v = loadColumns(2 * numAtoms + 2, item, "XRay_F1.txt", "anomalous Rayleigh data");
        _F2v = loadColumns(2 * numAtoms + 2, item, "XRay_F2.txt", "anomalous Rayleigh data");

        // precalculate the logarithms needed to interpolate the atomic form factors
        _logRatioFFv = logRatios(_FFv);

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

    ////////////////////////////////////////////////////////////////////

    double AnomalousRayleighHelper::sectionSca(double lambda, int Z, int /*N*/) const
    {
        // interpolate from table, and:
        // - below lower table limit: use Z^2 * Thomson scattering
        // - above upper table limit: use clamped value
        double x = scaledEnergy(lambda);
        if (x < _RSAv[2 * Z][0]) return Z * Z * Constants::sigmaThomson();
        return NR::clampedValue<NR::interpolateLogLog>(x, _RSAv[2 * Z], _RSAv[2 * Z + 1]);
    }

    ////////////////////////////////////////////////////////////////////

    double AnomalousRayleighHelper::phaseFunctionValue(double x, double costheta, int Z) const
    {
        constexpr double norm = 3. / 4. * Constants::sigmaThomson();
        double phase = 1. + costheta * costheta;
        double section = NR::clampedValue<NR::interpolateLogLog>(x, _RSAv[2 * Z], _RSAv[2 * Z + 1]);
        double sintheta2 = sqrt(0.5 * (1 - costheta));
        double form = interpolateQ(x, sintheta2, _FFv, _logRatioFFv, Z);
        double form1 = NR::clampedValue<NR::interpolateLogLin>(x, _F1v[2 * Z], _F1v[2 * Z + 1]);  // negative values
        double form2 = NR::clampedValue<NR::interpolateLogLog>(x, _F2v[2 * Z], _F2v[2 * Z + 1]);
        double formsum = form + form1;
        return norm / section * phase * (formsum * formsum + form2 * form2);
    }

    ////////////////////////////////////////////////////////////////////

    double AnomalousRayleighHelper::generateCosineFromPhaseFunction(double x, double Z) const
    {
        // construct the normalized cumulative phase function distribution for this x
        Array thetaXv;
        Array formv = interpolateQ(x, _sintheta2v, _FFv, _logRatioFFv, Z);
        NR::cdf(thetaXv, maxTheta, [this, x, Z, &formv](int t) {
            t += 1;
            double phase = 1. + _cos2thetav[t];
            double form = formv[t];
            double form1 = NR::clampedValue<NR::interpolateLogLin>(x, _F1v[2 * Z], _F1v[2 * Z + 1]);
            double form2 = NR::clampedValue<NR::interpolateLogLog>(x, _F2v[2 * Z], _F2v[2 * Z + 1]);
            double formsum = form + form1;
            return phase * (formsum * formsum + form2 * form2) * _sinthetav[t];
        });

        // draw a random cosine from this distribution
        return _random->cdfLinLin(_costhetav, thetaXv);
    }

    ////////////////////////////////////////////////////////////////////

    void AnomalousRayleighHelper::peeloffScattering(double& I, double& lambda, int Z, int /*N*/, Direction bfk,
                                                    Direction bfkobs) const
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

    ////////////////////////////////////////////////////////////////////

    Direction AnomalousRayleighHelper::performScattering(double& lambda, int Z, int /*N*/, Direction bfk) const
    {
        double x = scaledEnergy(lambda);

        // for low energies use Thomson scattering
        if (x < _RSAv[2 * Z][0]) return _dpf.performScattering(bfk, nullptr);

        // otherwise use Rayleigh scattering
        return _random->direction(bfk, generateCosineFromPhaseFunction(x, Z));
    }
}
