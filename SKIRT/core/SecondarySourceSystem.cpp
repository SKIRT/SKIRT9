/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SecondarySourceSystem.hpp"
#include "AngularDistributionInterface.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "Gas.hpp"
#include "Log.hpp"
#include "MaterialMix.hpp"
#include "MediumSystem.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PhotonPacket.hpp"
#include "PolarizationProfileInterface.hpp"
#include "ProbePhotonPacketInterface.hpp"
#include "ProcessManager.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "VelocityInterface.hpp"
#include "WavelengthDistribution.hpp"

////////////////////////////////////////////////////////////////////

SecondarySourceSystem::SecondarySourceSystem(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

void SecondarySourceSystem::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    _config = find<Configuration>(false);
    _ms = find<MediumSystem>(false);
    _random = find<Random>(false);
}

////////////////////////////////////////////////////////////////////

void SecondarySourceSystem::installLaunchCallBack(ProbePhotonPacketInterface* callback)
{
    if (_callback) throw FATALERROR("Cannot install more than one photon packet launch probe");
    _callback = callback;
}

////////////////////////////////////////////////////////////////////

// common actions for gas and dust, but they need to fill in different members
namespace
{
    // calculate the absorbed (and thus to be emitted) dust or gas luminosity for each spatial
    // cell this can be somewhat time-consuming, so we do this in parallel
    void calculateLv(Array& Lv, MaterialMix::MaterialType t, const MediumSystem* ms)
    {
        int numCells = ms->numCells();
        Lv.resize(numCells);
        ms->find<ParallelFactory>()->parallelDistributed()->call(numCells, [&](size_t firstIndex, size_t numIndices) {
            for (size_t m = firstIndex; m != firstIndex + numIndices; ++m) Lv[m] = ms->absorbedLuminosity(m, t);
        });
        ProcessManager::sumToAll(Lv);
    }

    // calculate the launch weight for each cell, given their luminosities and a composite bias fraction
    void calculateWv(Array& Wv, const Array& Lv, double xi, int numCells)
    {
        // determine a uniform weight for each cell with non-negligable emission, and normalize to unity
        Array wv(numCells);
        for (int m = 0; m != numCells; ++m) wv[m] = Lv[m] > 0 ? 1. : 0.;
        wv /= wv.sum();

        // calculate the final, composite-biased launch weight for each cell, normalized to unity
        Wv = (1 - xi) * Lv + xi * wv;
    }

    // calculate the binning to translate a packet history index to a (sorted) cell index p, given
    // the weight for each cell (Wv), the mapping from p to m (mv), and a number of photon packets
    // (Npp). An offset can be given in case both dust and gas need an Iv array.
    void calculateIv(vector<size_t>& Iv, const Array& Wv, const vector<int>& mv, size_t Npp, size_t offset,
                     int numCells)
    {
        Iv.resize(numCells + 1);
        Iv[0] = offset;
        double W = 0.;
        for (int p = 1; p != numCells; ++p)
        {
            // track the cumulative normalized weight as a floating point number and limit the
            // index to numPackets to avoid issues with rounding errors
            W += Wv[mv[p - 1]];
            Iv[p] = offset + min(Npp, static_cast<size_t>(std::round(W * Npp)));
        }
        Iv[numCells] = offset + Npp;
    }
}

////////////////////////////////////////////////////////////////////

bool SecondarySourceSystem::prepareForLaunch(size_t numPackets)
{
    // --------- luminosities 1 ---------
    calculateLv(_Ldustv, MaterialMix::MaterialType::Dust, _ms);
    calculateLv(_Lgasv, MaterialMix::MaterialType::Gas, _ms);

    // --------- library mapping ---------

    // obtain the spatial cell library mapping (from cell indices to library entry indices);
    // pass the cell luminosities to the library so it can avoid mapping zero-luminosity cells
    _nv = _config->cellLibrary()->mapping(_Ldustv);

    // construct a list of spatial cell indices sorted so that cells belonging to the same entry are consecutive
    int numCells = _ms->numCells();
    _mv.resize(numCells);
    for (int m = 0; m != numCells; ++m) _mv[m] = m;
    std::sort(begin(_mv), end(_mv), [this](int m1, int m2) { return _nv[m1] < _nv[m2]; });

    // --------- luminosities 2 ---------

    // force cells that have not been mapped by the library to zero luminosity
    // so that we won't launch photon packets for these cells
    for (int m = 0; m != numCells; ++m)
        if (_nv[m] < 0) _Ldustv[m] = 0.;

    // calculate the total luminosity; if both are zero, report failure
    _Ldust = _Ldustv.sum();
    _Lgas = _Lgasv.sum();

    // use these flags to deal with the special case that there is no dust or no gas
    _launchDust = _Ldust > 0;
    _launchGas = _Lgas > 0;

    // calculate the average luminosity contribution for each packet
    // TODO: implement option to do only dust or only gas, and variable ratio
    size_t Nppdust = 0;
    size_t Nppgas = 0;
    _Lppdust = 0.;
    _Lppgas = 0.;
    if (_launchDust && _launchGas)
    {
        Nppdust = numPackets / 2;
        Nppgas = numPackets - Nppdust;
        _Lppdust = _Ldust / Nppdust;
        _Lppgas = _Lgas / Nppgas;
    }
    else if (_launchDust)
    {
        Nppdust = numPackets;
        _Lppdust = _Ldust / Nppdust;
    }
    else if (_launchGas)
    {
        Nppgas = numPackets;
        _Lppgas = _Lgas / Nppgas;
    }
    else
    {
        return false;
    }

    // normalize the individual luminosities to unity
    _Ldustv /= _Ldust;
    _Lgasv /= _Lgas;

    // --------- weights ---------

    double xi = _config->secondarySpatialBias();

    // determine the first history index for each cell, using the adjusted cell ordering so that
    // all photon packets for a given library entry are launched consecutively
    if (_launchDust)
    {
        calculateWv(_Wdustv, _Ldustv, xi, numCells);
        calculateIv(_Idustv, _Wdustv, _mv, Nppdust, 0, numCells);
    }

    // if there were 0 dust packets, start from 0. if there was 1, start from index 1, etc. If not
    // _launchGas, then this index == Nppdust == numPackets. The history index for the packets only
    // goes up to numPackets, so no gas packets will be launched, as intended.
    _startGasIndex = Nppdust;

    if (_launchGas)
    {
        calculateWv(_Wgasv, _Lgasv, xi, numCells);
        calculateIv(_Igasv, _Wgasv, _mv, Nppgas, _startGasIndex, numCells);
    }

    // --------- logging ---------

    auto log = find<Log>();

    // spatial cells
    int emittingCells = 0;  // number of nonzero luminosity cells
    for (int m = 0; m != numCells; ++m)
        if (_Ldustv[m] > 0. || _Lgasv[m] > 0.) emittingCells++;
    log->info("Emitting from " + std::to_string(emittingCells) + " out of " + std::to_string(numCells)
              + " spatial cells");

    // library entries
    int numEntries = _config->cellLibrary()->numEntries();
    vector<int> mapped(numEntries);  // number of cells mapped to each library entry
    for (int m = 0; m != numCells; ++m)
        if (_nv[m] >= 0) mapped[_nv[m]]++;

    int usedEntries = 0;     // number of library entries that have at least one mapped cell
    int maxMappedCells = 0;  // largest number of cells mapped to a library entry
    int totMappedCells = 0;  // total number of cells mapped to a library entry
    for (int n = 0; n != numEntries; ++n)
    {
        if (mapped[n] > 0)
        {
            usedEntries++;
            maxMappedCells = max(mapped[n], maxMappedCells);
            totMappedCells += mapped[n];
        }
    }

    if (numEntries != numCells || totMappedCells != numCells || numEntries != usedEntries || maxMappedCells != 1)
    {
        log->info("  Using " + std::to_string(usedEntries) + " out of " + std::to_string(numEntries)
                  + " library entries");
        log->info("  Largest number of cells per library entry: " + std::to_string(maxMappedCells));
        log->info("  Average number of cells per (used) library entry: "
                  + StringUtils::toString(static_cast<double>(totMappedCells) / usedEntries, 'f', 1));
    }

    // report success
    return true;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // An instance of this class obtains and/or calculates the information needed to launch photon packets
    // from the dust in a given cell in the spatial grid, and remembers the information for fast retrieval.
    // This information includes the normalized regular and cumulative dust emission spectra, calculated from
    // the stored radiation field and the dust properties using the configured emission calculator, and
    // the average bulk velocity of the material in the cell, obtained from the medium system.
    class DustCellEmission : public VelocityInterface
    {
    private:
        // information initialized once, during the first call to calculateIfNeeded()
        MediumSystem* _ms{nullptr};  // the medium system
        Array _wavelengthGrid;       // the dust emission wavelength grid
        Range _wavelengthRange;      // the range of the dust emission wavelength grid
        int _numWavelengths{0};      // the number of wavelengths in the dust emission wavelength grid
        vector<int> _hv;             // a list of the media indices for the media containing dust
        int _numMedia{0};            // the number of dust media in the system (and thus the size of hv)
        int _numCells{0};            // the number of cells in the spatial grid (and thus the size of mv and nv)

        // information on a particular spatial cell, initialized by calculateIfNeeded()
        int _p{-1};                // spatial cell launch-order index
        int _n{-1};                // library entry index
        vector<Array> _evv;        // emissivity spectrum for each medium component, if applicable
        Array _lambdav, _pv, _Pv;  // normalized emission spectrum
        Vec _bfv;                  // bulk velocity

    public:
        DustCellEmission() {}

        // calculates the emission information for the given cell if it is different from what's already stored
        //   p:  launch-order cell index (cells mapped to a given library entry have consecutive p indices)
        //   mv: map from launch-order cell index p to regular cell index m
        //   nv: map from regular cell index m to library entry index n
        //   ms: medium system
        //   config: configuration object
        void calculateIfNeeded(int p, const vector<int>& mv, const vector<int>& nv, MediumSystem* ms,
                               Configuration* config)
        {
            // if this photon packet is launched from the same cell as the previous one, we don't need to do anything
            if (p == _p) return;

            // when called for the first time, construct a list of dust media and cache some other info
            if (_p == -1)
            {
                _ms = ms;
                auto wavelengthGrid = config->dustEmissionWLG();
                _wavelengthGrid = wavelengthGrid->extlambdav();
                _wavelengthRange = wavelengthGrid->wavelengthRange();
                _numWavelengths = _wavelengthGrid.size();
                for (int h = 0; h != ms->numMedia(); ++h)
                    if (ms->isDust(h)) _hv.push_back(h);
                _numMedia = _hv.size();
                _numCells = ms->numCells();
                _evv.resize(ms->numMedia());
            }

            // remember the new cell index and map to the other indices
            _p = p;
            int m = mv[p];
            int n = nv[m];

            // if this new cell maps to a new library entry, we need to process the library entry
            if (n != _n)
            {
                // remember the new library entry index
                _n = n;

                // determine the number of cells mapped to this library entry (they are consecutive in p)
                int pp = p + 1;
                for (; pp != _numCells; ++pp)
                    if (nv[mv[pp]] != n) break;
                int numMappedCells = pp - p;

                // if only a single cell maps to the library entry, we can simply calculate its emission
                if (numMappedCells == 1)
                {
                    calculateSingleSpectrum(_ms->meanIntensity(m), m);
                }

                // if multiple cells map to the library entry, we use the average radiation field for these cells
                else
                {
                    Array Jv = _ms->meanIntensity(m);
                    for (int i = 1; i != numMappedCells; ++i) Jv += _ms->meanIntensity(mv[p + i]);
                    Jv /= numMappedCells;

                    // if there is a single dust medium (and assuming that there are no variable dust mixes),
                    // we can use a single emission spectrum for all cells mapped to the library entry, calculated
                    // using the average radiation field, because the cells differ only in dust density, which is
                    // irrelevant because the emission spectrum is normalized anyway
                    if (_numMedia == 1)
                    {
                        calculateSingleSpectrum(Jv, m);
                    }

                    // otherwise, we need to calculate and remember the emission spectrum for each medium component
                    // (still using the average radiation field and assuming that there are no variable dust mixes)
                    // so that we can apply the relative density weights for each cell later on
                    else
                    {
                        calculateEmissivityPerMedium(Jv, m);
                        calculateWeightedSpectrum(m);
                    }
                }
            }

            // if this new cell maps to the same library entry, we can use its cached info
            else
            {
                // if there is a single dust medium we can use the already calculated emission spectrum;
                // otherwise, we need to apply the relative density weights for this cell to the calculated
                // emission spectra for each dust medium, and renormalize the resulting spectrum
                if (_numMedia != 1)
                {
                    calculateWeightedSpectrum(m);
                }
            }

            // remember the average bulk velocity for this cell
            _bfv = ms->bulkVelocity(m);
        }

    private:
        // calculate the emission spectrum for the specified radiation field and the dust mixes of the specified cell,
        // and store the result in the data members _lambdav, _pv, _Pv
        void calculateSingleSpectrum(const Array& Jv, int m)
        {
            // accumulate the emmissivity spectrum for all dust medium components in the cell, weighted by density
            Array ev(_numWavelengths);
            for (int h : _hv) ev += _ms->numberDensity(m, h) * _ms->mix(m, h)->emissivity(Jv);

            // calculate the normalized plain and cumulative distributions
            NR::cdf<NR::interpolateLogLog>(_lambdav, _pv, _Pv, _wavelengthGrid, ev, _wavelengthRange);
        }

        // calculate the emmissivity spectra for the specified radiation field and the dust mixes of the specified cell,
        // and store the individual spectra in the data members _evv
        void calculateEmissivityPerMedium(const Array& Jv, int m)
        {
            for (int h : _hv) _evv[h] = _ms->mix(m, h)->emissivity(Jv);
        }

        // calculate the emission spectrum for the specified cell, weighted across multiple media by density,
        // given the precalculated emissivity spectra _evv for each medium,
        // and store the result in the data members _lambdav, _pv, _Pv
        void calculateWeightedSpectrum(int m)
        {
            // accumulate the emmissivity spectrum for all dust medium components in the cell, weighed by density
            Array ev(_numWavelengths);
            for (int h : _hv) ev += _ms->numberDensity(m, h) * _evv[h];

            // calculate the normalized plain and cumulative distributions
            NR::cdf<NR::interpolateLogLog>(_lambdav, _pv, _Pv, _wavelengthGrid, ev, _wavelengthRange);
        }

    public:
        // returns a random wavelength generated from the spectral distribution
        double generateWavelength(Random* random) const { return random->cdfLogLog(_lambdav, _pv, _Pv); }

        // returns the normalized specific luminosity for the given wavelength
        double specificLuminosity(double lambda) const
        {
            return NR::value<NR::interpolateLogLog>(lambda, _lambdav, _pv);
        }

        Vec velocity() const override { return _bfv; }
    };

    // setup a DustCellEmission instance for each parallel execution thread to cache dust emission information
    thread_local DustCellEmission t_dustcell;

    // Same as the above, but for gas. Some things are similar, others are not. Implement with some
    // duplication for now.
    class GasCellEmission : public VelocityInterface
    {
    private:
        MediumSystem* _ms{nullptr};
        Array _wavelengthGrid;
        Range _wavelengthRange;

        int _p{-1};
        Array _lambdav, _pv, _Pv;
        Vec _bfv;

    public:
        GasCellEmission() {}

        void calculateIfNeeded(int p, const vector<int>& mv, MediumSystem* ms)
        {
            if (p == _p) return;
            if (_p == -1)
            {
                _ms = ms;
                _wavelengthGrid = Gas::emissivityLambdav();
                _wavelengthRange = Range(_wavelengthGrid[0], _wavelengthGrid[_wavelengthGrid.size() - 1]);
            }
            _p = p;
            int m = mv[p];
            calculateSingleSpectrum(m);
            _bfv = ms->bulkVelocity(m);
        }

    private:
        void calculateSingleSpectrum(int m)
        {
            Array ev = Gas::emissivity(m);
            NR::cdf<NR::interpolateLogLog>(_lambdav, _pv, _Pv, _wavelengthGrid, ev, _wavelengthRange);
        }

    public:
        double generateWavelength(Random* random) const { return random->cdfLogLog(_lambdav, _pv, _Pv); }
        double specificLuminosity(double lambda) const
        {
            return NR::value<NR::interpolateLogLog>(lambda, _lambdav, _pv);
        }
        Vec velocity() const override { return _bfv; }
    };

    thread_local GasCellEmission t_gascell;

    // An instance of this class obtains the information needed to determine
    // angular distribution probabilities and polarization components for
    // photon packets emitted from the dust in a given cell in the spatial grid,
    // and remembers some information for fast retrieval.
    class DustCellPolarisedEmission : public AngularDistributionInterface, public PolarizationProfileInterface
    {
    private:
        // information initialized once, during the first call to calculateIfNeeded()
        MediumSystem* _ms{nullptr};  // the medium system
        vector<int> _hv;             // a list of the media indices for the media containing dust
        int _numMedia{0};            // the number of dust media in the system (and thus the size of hv)

        // information on a particular spatial cell, initialized by calculateIfNeeded()
        int _m{-1};        // spatial cell index
        Vec _B_direction;  // direction of the magnetic field

        // information on a particular photon packet, initialized by calculate()
        double _lambda{-1.};  // wavelength of the last photon packet

    public:
        DustCellPolarisedEmission() {}

        // calculates the emission information for the given cell if it is different from what's already stored
        //   m:  cell index
        //   ms: medium system
        //   config: configuration object
        void calculateIfNeeded(int m, MediumSystem* ms, Configuration* config)
        {
            // if this photon packet is launched from the same cell as the previous one, we don't need to do anything
            if (!config->hasSpheroidalPolarization() || m == _m) return;

            // when called for the first time, construct a list of dust media and cache some other info
            if (_m == -1)
            {
                _ms = ms;
                for (int h = 0; h != ms->numMedia(); ++h)
                    if (ms->isDust(h)) _hv.push_back(h);
                _numMedia = _hv.size();
            }

            // remember the new cell index and map to the other indices
            _m = m;

            // normalise the magnetic field direction
            _B_direction = _ms->magneticField(_m);
            _B_direction /= _B_direction.norm();
        }

        // calculate the emission information that is always out of date
        // in practice, we use this function to set the wavelength of the emission
        //   lambda: wavelength
        void calculate(double lambda) { _lambda = lambda; }

        // this AngularDistributionInterface implementation returns the probability for the given direction
        //   bfk: direction
        double probabilityForDirection(Direction bfk) const override
        {
            // get the zenith angle of the direction
            // this is the angle between the magnetic field direction (z axis of the reference frame)
            // and the direction
            const double theta = std::acos(Vec::dot(_B_direction, bfk));

            // compute
            //   - Qabstheta: the absorption coefficient for the direction
            //   - Qabstot: the agnular averaged absorption coefficient
            double Qabstheta = 0.;
            double Qabstot = 0.;
            // we need to sum over all dust species
            for (int h : _hv)
            {
                const Array& thetas = _ms->mix(_m, h)->thetaGrid();
                const Array& Qabs = _ms->mix(_m, h)->sectionsAbs(_lambda);
                Qabstheta += _ms->numberDensity(_m, h) * NR::value<NR::interpolateLinLin>(theta, thetas, Qabs);
                Qabstot += _ms->numberDensity(_m, h) * _ms->mix(_m, h)->sectionAbs(_lambda);
            }
            // the probability for this direction is the ratio of the two values
            return Qabstheta / Qabstot;
        }

        // generate a random direction based on the anisotropic emission profile
        Direction generateDirection(Random* random) const
        {
            // first generate a random azimuth angle
            const double phi = 2. * M_PI * random->uniform();
            // compute its sine and cosine
            const double cosphi = cos(phi);
            const double sinphi = sin(phi);

            // now generate a random zenith angle
            // first, compute the cumulative distribution for a photon at the given wavelength
            // create an empty array with the same size as the zenith angle grid
            const Array& thetas = _ms->mix(_m, _hv[0])->thetaGrid();
            Array cdf(thetas.size());
            // add the distributions for the different components
            for (int h : _hv)
            {
                const Array& Qabs = _ms->mix(_m, h)->sectionsAbs(_lambda);
                cdf += _ms->numberDensity(_m, h) * Qabs;
            }
            // convert to a cumulative distribution
            for (size_t i = 1; i < thetas.size(); ++i)
            {
                cdf[i] += cdf[i - 1];
            }
            // normalise the distribution
            for (size_t i = 0; i < thetas.size(); ++i)
            {
                cdf[i] /= cdf[thetas.size() - 1];
            }
            // draw a random uniform deviate
            // get the corresponding zenith angle
            const double theta = random->cdfLinLin(thetas, cdf);
            // compute the sine and cosine
            const double costheta = cos(theta);
            const double sintheta = sin(theta);

            // generate a random direction
            return Direction(sintheta * cosphi, sintheta * sinphi, costheta);
        }

        // this PolarizationProfileInterface implementation returns the Stokes vector for the given emission direction
        //   bfk: emission direction
        StokesVector polarizationForDirection(Direction bfk) const override
        {
            // get the zenith angle of the direction
            // this is the angle between the magnetic field direction (z axis of the reference frame)
            // and the direction
            const double theta = std::acos(Vec::dot(_B_direction, bfk));

            // we sum over all dust species
            double Qabsval = 0.;
            double Qabspolval = 0.;
            for (int h : _hv)
            {
                const Array& thetas = _ms->mix(_m, h)->thetaGrid();
                const Array& Qabs = _ms->mix(_m, h)->sectionsAbs(_lambda);
                const Array& Qabspol = _ms->mix(_m, h)->sectionsAbspol(_lambda);
                Qabsval += _ms->numberDensity(_m, h) * NR::value<NR::interpolateLinLin>(theta, thetas, Qabs);
                Qabspolval += _ms->numberDensity(_m, h) * NR::value<NR::interpolateLinLin>(theta, thetas, Qabspol);
            }
            // the reference direction is any direction perpendicular to the
            // propagation direction and the magnetic field direction
            return StokesVector(Qabsval, Qabspolval, 0., 0., Direction(Vec::cross(bfk, _B_direction)));
        }
    };

    // setup a DustCellPolarisedEmission instance for each parallel execution thread to cache dust emission information
    thread_local DustCellPolarisedEmission t_dustcellpol;
}

////////////////////////////////////////////////////////////////////

void SecondarySourceSystem::launch(PhotonPacket* pp, size_t historyIndex) const
{
    // select the spatial cell from which to launch based on the history index of this photon packet
    bool isDust = historyIndex < _startGasIndex;

    // TODO: remove this check once everything works
    if ((isDust && !_launchDust) || (!isDust && !_launchGas)) FATALERROR("Something went wrong in prepareForLaunch");

    size_t p = isDust ? std::upper_bound(_Idustv.cbegin(), _Idustv.cend(), historyIndex) - _Idustv.cbegin() - 1
                      : std::upper_bound(_Igasv.cbegin(), _Igasv.cend(), historyIndex) - _Igasv.cbegin() - 1;

    auto m = _mv[p];

    if (isDust)
    {
        // calculate the emission spectrum and bulk velocity for this cell, if not already available
        t_dustcell.calculateIfNeeded(p, _mv, _nv, _ms, _config);
        t_dustcellpol.calculateIfNeeded(m, _ms, _config);
    }
    else
    {
        t_gascell.calculateIfNeeded(p, _mv, _ms);
    }

    // generate a random wavelength from the emission spectrum for the cell and/or from the bias distribution
    double lambda, w;
    {
        double xi = _config->secondaryWavelengthBias();
        if (!xi)
        {
            // no biasing -- simply use the intrinsic spectral distribution
            lambda = isDust ? t_dustcell.generateWavelength(_random) : t_gascell.generateWavelength(_random);
            w = 1.;
        }
        else
        {
            // biasing -- use one or the other distribution
            if (_random->uniform() > xi)
                lambda = isDust ? t_dustcell.generateWavelength(_random) : t_gascell.generateWavelength(_random);
            else
                lambda = _config->secondaryWavelengthBiasDistribution()->generateWavelength();

            // calculate the compensating weight factor
            double s = isDust ? t_dustcell.specificLuminosity(lambda) : t_gascell.specificLuminosity(lambda);
            if (!s)
            {
                // if the wavelength can't occur in the intrinsic distribution,
                // the weight factor is zero regardless of the probability in the bias distribution
                // (handling this separately also avoids NaNs in the pathological case s=b=0)
                w = 0.;
            }
            else
            {
                // regular composite bias weight
                double b = _config->secondaryWavelengthBiasDistribution()->probability(lambda);
                w = s / ((1 - xi) * s + xi * b);
            }
        }
    }

    // get the weighted luminosity corresponding to this cell
    double L = isDust ? _Lppdust * _Ldustv[m] / _Wdustv[m] : _Lppgas * _Lgasv[m] / _Wgasv[m];

    // generate a random position in this spatial cell
    Position bfr = _ms->grid()->randomPositionInCell(m);

    // provide a redshift interface for the appropriate velocity, if it is nonzero
    VelocityInterface* bvi = nullptr;
    if (isDust && !t_dustcell.velocity().isNull())
        bvi = &t_dustcell;
    else if (!t_gascell.velocity().isNull())
        bvi = &t_gascell;

    // provide a polarisation interface for polarised emission, if applicable
    // generate a random emission direction
    DustCellPolarisedEmission* dpe = nullptr;
    Direction bfk;
    if (isDust && _config->hasSpheroidalPolarization())
    {
        t_dustcellpol.calculate(lambda);
        dpe = &t_dustcellpol;
        bfk = t_dustcellpol.generateDirection(_random);
    }
    else
    {
        bfk = _random->direction();
    }

    pp->launch(historyIndex, lambda, L * w, bfr, bfk, bvi, dpe, dpe);

    // add origin info (we combine all medium components, so we cannot differentiate them here)
    pp->setSecondaryOrigin(0);

    // invoke launch call-back if installed
    if (_callback) _callback->probePhotonPacket(pp);
}

////////////////////////////////////////////////////////////////////
