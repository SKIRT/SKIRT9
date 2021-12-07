/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustSecondarySource.hpp"
#include "AngularDistributionInterface.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PhotonPacket.hpp"
#include "PolarizationProfileInterface.hpp"
#include "ProcessManager.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "Units.hpp"
#include "VelocityInterface.hpp"
#include "WavelengthDistribution.hpp"

////////////////////////////////////////////////////////////////////

double DustSecondarySource::prepareLuminosities()
{
    // cache some pointers for later use
    _config = find<Configuration>();
    _ms = find<MediumSystem>();
    _random = find<Random>();

    int numCells = _ms->numCells();

    // --------- luminosities 1 ---------

    // calculate the absorbed (and thus to be emitted) dust luminosity for each spatial cell
    // this can be somewhat time-consuming, so we do this in parallel
    _Lv.resize(numCells);
    find<ParallelFactory>()->parallelDistributed()->call(numCells, [this](size_t firstIndex, size_t numIndices) {
        for (size_t m = firstIndex; m != firstIndex + numIndices; ++m)
        {
            _Lv[m] = _ms->dustLuminosity(m);
        }
    });
    ProcessManager::sumToAll(_Lv);

    // --------- library mapping ---------

    // obtain the spatial cell library mapping (from cell indices to library entry indices);
    // pass the cell luminosities to the library so it can avoid mapping zero-luminosity cells
    _nv = _config->cellLibrary()->mapping(_Lv);

    // construct a list of spatial cell indices sorted so that cells belonging to the same entry are consecutive
    _mv.resize(numCells);
    for (int m = 0; m != numCells; ++m) _mv[m] = m;
    std::sort(begin(_mv), end(_mv), [this](int m1, int m2) { return _nv[m1] < _nv[m2]; });

    // --------- luminosities 2 ---------

    // force cells that have not been mapped by the library to zero luminosity
    // so that we won't launch photon packets for these cells
    for (int m = 0; m != numCells; ++m)
        if (_nv[m] < 0) _Lv[m] = 0.;

    // calculate  the total luminosity, and normalize the individual luminosities to unity
    double L = _Lv.sum();
    _Lv /= L;

    // --------- logging ---------

    auto log = find<Log>();
    auto units = find<Units>();

    // luminosity and spatial cells
    int emittingCells = 0;  // number of nonzero luminosity cells
    for (int m = 0; m != numCells; ++m)
        if (_Lv[m] > 0.) emittingCells++;
    log->info("Dust luminosity: " + StringUtils::toString(units->obolluminosity(L), 'g') + " "
              + units->ubolluminosity());
    log->info("  Emitting from " + std::to_string(emittingCells) + " out of " + std::to_string(numCells)
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

    // return the total luminosity
    return L;
}

////////////////////////////////////////////////////////////////////

void DustSecondarySource::preparePacketMap(size_t firstIndex, size_t numIndices)
{
    int numCells = _ms->numCells();

    // --------- weights ---------

    // determine a uniform weight for each cell with non-negligable emission, and normalize to unity
    Array wv(numCells);
    for (int m = 0; m != numCells; ++m) wv[m] = _Lv[m] > 0 ? 1. : 0.;
    wv /= wv.sum();

    // calculate the final, composite-biased launch weight for each cell, normalized to unity
    double xi = _config->secondarySpatialBias();
    _Wv = (1 - xi) * _Lv + xi * wv;

    // determine the first history index for each cell, using the adjusted cell ordering so that
    // all photon packets for a given library entry are launched consecutively
    _Iv.resize(numCells + 1);
    _Iv[0] = firstIndex;
    double W = 0.;
    for (int p = 1; p != numCells; ++p)
    {
        // track the cumulative normalized weight as a floating point number
        // and limit the index to firstIndex+numIndices to avoid issues with rounding errors
        W += _Wv[_mv[p - 1]];
        _Iv[p] = firstIndex + min(numIndices, static_cast<size_t>(std::round(W * numIndices)));
    }
    _Iv[numCells] = firstIndex + numIndices;
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
        // instances of this class are allocated in thread-local storage, which means that the constructor is
        // invoked under the hood for each new execution thread; hence it does trivial initialization only
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
            // when called for the first time for a given simulation, cache some info
            if (_ms != ms)
            {
                _p = -1;
                _n = -1;
                _ms = ms;
                auto wavelengthGrid = config->dustEmissionWLG();
                _wavelengthGrid = wavelengthGrid->extlambdav();
                _wavelengthRange = wavelengthGrid->wavelengthRange();
                _numWavelengths = _wavelengthGrid.size();
                _hv = ms->dustMediumIndices();
                _numMedia = _hv.size();
                _numCells = ms->numCells();
                _evv.resize(ms->numMedia());
            }

            // if this photon packet is launched from the same cell as the previous one, we don't need to do anything
            if (p == _p) return;

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
                    calculateSingleSpectrum(m);
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
        // calculate the emission spectrum for the dust mixes of the specified cell,
        // and store the result in the data members _lambdav, _pv, _Pv
        void calculateSingleSpectrum(int m)
        {
            // get the emmissivity spectrum for all dust medium components in the cell
            const Array& ev = _ms->dustEmissionSpectrum(m);

            // calculate the normalized plain and cumulative distributions
            NR::cdf<NR::interpolateLogLog>(_lambdav, _pv, _Pv, _wavelengthGrid, ev, _wavelengthRange);
        }

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

    // An instance of this class obtains the information needed to determine angular distribution probabilities
    // and polarization components for photon packets emitted from the dust in a given cell in the spatial grid,
    // and remembers some information for fast retrieval. It is actually used only when there are spheroidal
    // dust grains in the simulation.
    class DustCellPolarisedEmission : public AngularDistributionInterface, public PolarizationProfileInterface
    {
    private:
        // information initialized once, during the first call to calculateIfNeeded()
        MediumSystem* _ms{nullptr};  // the medium system
        vector<int> _hv;             // a list of the media indices for the media containing dust

        // information on a particular spatial cell, initialized by calculateIfNeeded()
        int _m{-1};        // spatial cell index
        Vec _B_direction;  // direction of the magnetic field

        // information on a particular photon packet, initialized by calculateIfNeeded()
        double _lambda{-1.};  // wavelength of the photon packet

    public:
        // instances of this class are allocated in thread-local storage, which means that the constructor is
        // invoked under the hood for each new execution thread; hence it does trivial initialization only
        DustCellPolarisedEmission() {}

        // calculates the emission information for the given cell and photon packet wavelength
        // if it is different from what's already stored
        void calculateIfNeeded(int m, MediumSystem* ms, double lambda)
        {
            // always remember the photon packet wavelength
            _lambda = lambda;

            // if this packet is launched from the same cell as the previous one, we don't need to do anything else
            if (m == _m) return;

            // when called for the first time, cache some info
            if (_m == -1)
            {
                _ms = ms;
                _hv = ms->dustMediumIndices();
            }

            // remember the new cell index and map to the other indices
            _m = m;

            // normalise the magnetic field direction
            _B_direction = _ms->magneticField(_m);
            _B_direction /= _B_direction.norm();
        }

        // this AngularDistributionInterface implementation returns the probability for the given direction
        double probabilityForDirection(Direction bfk) const override
        {
            // get the zenith angle of the direction
            // this is the angle between the magnetic field direction (z axis of the reference frame)
            // and the direction
            const double theta = std::acos(Vec::dot(_B_direction, bfk));

            // compute
            //   - Qabstheta: the absorption coefficient for the direction
            //   - Qabstot: the angular averaged absorption coefficient
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

            // now generate a random zenith angle
            // first, compute the cumulative distribution for a photon at the given wavelength
            // create empty PDF and CDF arrays with the same size as the zenith angle grid
            const Array& thetas = _ms->mix(_m, _hv[0])->thetaGrid();
            Array pdf(thetas.size()), cdf(thetas.size());
            // add the distributions for the different components
            for (int h : _hv)
            {
                const Array& Qabs = _ms->mix(_m, h)->sectionsAbs(_lambda);
                pdf += _ms->numberDensity(_m, h) * Qabs;
            }
            // the CDF is computed from the PDF by integration
            // since this is an integration is spherical coordinates, we need to
            // multiply the PDF with the volume element sin(theta)
            for (size_t i = 0; i < thetas.size(); ++i)
            {
                pdf[i] *= std::sin(thetas[i]);
            }
            // convert to the CDF
            // note that we use the midpoint of each bin as integration point
            // and omit the common factors from the integration
            for (size_t i = 1; i < thetas.size(); ++i)
            {
                cdf[i] = cdf[i - 1] + pdf[i - 1] + pdf[i];
            }
            // normalise the CDF
            for (size_t i = 0; i < thetas.size(); ++i)
            {
                cdf[i] /= cdf[thetas.size() - 1];
            }
            // draw a random zenith angle from the CDF
            const double theta = random->cdfLinLin(thetas, cdf);

            // generate a random direction
            const Direction krand(theta, phi);

            // The random direction is defined in a frame that has the B direction
            // as z-axis (the x and y axis can be chosen arbitrarily because of the
            // symmetry in phi)
            // We need to (passively) rotate the random direction to the lab frame
            // where the z-axis is the actual z-axis
            // This corresponds to a rotation over -beta degrees along an axis
            // that is perpendicular to the plane through the B direction and the
            // z-axis, where beta is the angle between the B direction and the
            // z-axis
            // Note that in this case, it does not matter whether we rotate over
            // beta or -beta degrees, since the difference between these two is
            // simply a rotation over 180 degrees in phi, and we still sample a
            // uniform distribution in phi.
            // If the angle is small enough, we simply assume that B already lies
            // along the z-axis and return immediately.
            const double cosbeta = _B_direction.z();
            if (fabs(cosbeta) < 0.99999)
            {
                // B is not the z-axis. Compute beta and sin(beta)
                const double beta = std::acos(cosbeta);
                const double sinbeta = std::sin(beta);
                // compute the direction of the rotation axis, z x B
                const Direction n(Vec::cross(Vec(0., 0., 1.), _B_direction));
                // now use Rodrigues' rotation formula to carry out the rotation
                const Direction krand2(krand * cosbeta + Vec::cross(n, krand) * sinbeta
                                       + n * Vec::dot(n, krand) * (1. - cosbeta));
                return krand2;
            }
            else
            {
                return krand;
            }
        }

        // this PolarizationProfileInterface implementation returns the Stokes vector for the given emission direction
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

    // setup instances of the above classes to cache dust emission information for each parallel execution thread
    thread_local DustCellEmission t_dustcell;
    thread_local DustCellPolarisedEmission t_dustcellpol;
}

////////////////////////////////////////////////////////////////////

void DustSecondarySource::launch(PhotonPacket* pp, size_t historyIndex, double L) const
{
    // select the spatial cell from which to launch based on the history index of this photon packet
    auto p = std::upper_bound(_Iv.cbegin(), _Iv.cend(), historyIndex) - _Iv.cbegin() - 1;
    auto m = _mv[p];

    // calculate the weight related to biased source selection
    double ws = _Lv[m] / _Wv[m];

    // calculate the emission spectrum and bulk velocity for this cell, if not already available
    t_dustcell.calculateIfNeeded(p, _mv, _nv, _ms, _config);

    // generate a random wavelength from the emission spectrum for the cell and/or from the bias distribution
    double lambda, w;
    {
        double xi = _config->dustEmissionWavelengthBias();
        if (!xi)
        {
            // no biasing -- simply use the intrinsic spectral distribution
            lambda = t_dustcell.generateWavelength(_random);
            w = 1.;
        }
        else
        {
            // biasing -- use one or the other distribution
            if (_random->uniform() > xi)
                lambda = t_dustcell.generateWavelength(_random);
            else
                lambda = _config->dustEmissionWavelengthBiasDistribution()->generateWavelength();

            // calculate the compensating weight factor
            double s = t_dustcell.specificLuminosity(lambda);
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
                double b = _config->dustEmissionWavelengthBiasDistribution()->probability(lambda);
                w = s / ((1 - xi) * s + xi * b);
            }
        }
    }

    // generate a random position in this spatial cell
    Position bfr = _ms->grid()->randomPositionInCell(m);

    // provide a redshift interface for the appropriate velocity, if it is nonzero
    VelocityInterface* bvi = t_dustcell.velocity().isNull() ? nullptr : &t_dustcell;

    // provide a polarisation interface for polarised emission, if applicable
    // generate a random emission direction
    DustCellPolarisedEmission* dpe = nullptr;
    Direction bfk;
    if (_config->hasSpheroidalPolarization())
    {
        t_dustcellpol.calculateIfNeeded(m, _ms, lambda);
        dpe = &t_dustcellpol;
        bfk = t_dustcellpol.generateDirection(_random);
    }
    else
    {
        bfk = _random->direction();
    }

    pp->launch(historyIndex, lambda, L * ws * w, bfr, bfk, bvi, dpe, dpe);
}

////////////////////////////////////////////////////////////////////
