/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedSource.hpp"
#include "Band.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PhotonPacket.hpp"
#include "ProcessManager.hpp"
#include "Random.hpp"
#include "SEDFamily.hpp"
#include "Snapshot.hpp"
#include "VelocityInterface.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // maximum number of luminosity calculations between two invocations of infoIfElapsed()
    const size_t logProgressChunkSize = 10000;
}

////////////////////////////////////////////////////////////////////

void ImportedSource::setupSelfAfter()
{
    Source::setupSelfAfter();

    auto config = find<Configuration>();
    _oligochromatic = config->oligochromatic();

    // warn the user if this source's intrinsic wavelength range does not fully cover the configured wavelength range
    informAvailableWavelengthRange(_sedFamily->intrinsicWavelengthRange(), _sedFamily->type());

    // determine the wavelength range for this soource
    _wavelengthRange = config->sourceWavelengthRange();
    _wavelengthRange.intersect(_sedFamily->intrinsicWavelengthRange());
    if (_wavelengthRange.empty())
        throw FATALERROR("Intrinsic SED family wavelength range does not overlap source wavelength range");

    // cache other wavelength information depending on whether this is an oligo- or panchromatic simulation
    if (_oligochromatic)
    {
        _arbitaryWavelength = config->wavelengthGrid(nullptr)->wavelength(0);
        _xi = 1.;  // always use bias distribution for oligochromatic simulations
        _biasDistribution = config->oligoWavelengthBiasDistribution();
    }
    else
    {
        _arbitaryWavelength = _wavelengthRange.mid();
        _xi = wavelengthBias();
        _biasDistribution = wavelengthBiasDistribution();
    }

    // create the snapshot with preconfigured spatial columns
    _snapshot = createAndOpenSnapshot();

    // add optional columns if applicable
    if (!_oligochromatic && _importVelocity)
    {
        _snapshot->importVelocity();
        if (_importVelocityDispersion) _snapshot->importVelocityDispersion();
    }
    if (_importCurrentMass) _snapshot->importCurrentMass();
    _snapshot->importParameters(_sedFamily->parameterInfo());

    // notify about building search data structures if needed
    if (find<Configuration>()->snapshotsNeedGetEntities()) _snapshot->setNeedGetEntities();

    // read the data from file
    _snapshot->readAndClose();

    // construct a vector with the (normalized) luminosity for each entity
    int M = _snapshot->numEntities();
    if (M)
    {
        // integrating over the SED for each entity can be time-consuming, so we do this in parallel
        _Lv.resize(M);
        auto log = find<Log>();
        log->info("Calculating luminosities for " + std::to_string(M) + " imported entities...");
        log->infoSetElapsed(M);
        find<ParallelFactory>()->parallelDistributed()->call(M, [this, log](size_t firstIndex, size_t numIndices) {
            Array lambdav, pv, Pv;  // the contents of these arrays is not used, so this could be optimized if needed
            Array params;

            while (numIndices)
            {
                size_t currentChunkSize = min(logProgressChunkSize, numIndices);
                for (size_t m = firstIndex; m != firstIndex + currentChunkSize; ++m)
                {
                    _snapshot->parameters(m, params);
                    _Lv[m] = _sedFamily->cdf(lambdav, pv, Pv, _wavelengthRange, params);
                }
                log->infoIfElapsed("Calculated luminosities: ", currentChunkSize);
                firstIndex += currentChunkSize;
                numIndices -= currentChunkSize;
            }
        });
        ProcessManager::sumToAll(_Lv);

        // remember the total luminosity and normalize the vector
        _L = _Lv.sum();
        if (_L) _Lv /= _L;
    }
}

////////////////////////////////////////////////////////////////////

ImportedSource::~ImportedSource()
{
    delete _snapshot;
}

////////////////////////////////////////////////////////////////////

int ImportedSource::dimension() const
{
    return 3;
}

////////////////////////////////////////////////////////////////////

bool ImportedSource::hasVelocity() const
{
    if (_importVelocity)
    {
        // refuse velocity for oligochromatic simulations
        // (this function is called from Configure so we cannot precompute this during setup)
        auto config = find<Configuration>();
        if (!config->oligochromatic()) return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////

Range ImportedSource::wavelengthRange() const
{
    // don't rely on the cached _wavelengthRange because this function may be called during setup
    // *before* setupSelfAfter() has been executed (which initializes _wavelengthRange)
    auto config = find<Configuration>();
    auto wavelengthRange = config->sourceWavelengthRange();
    wavelengthRange.intersect(_sedFamily->intrinsicWavelengthRange());
    return wavelengthRange;
}

////////////////////////////////////////////////////////////////////

double ImportedSource::luminosity() const
{
    return _L;
}

////////////////////////////////////////////////////////////////////

double ImportedSource::specificLuminosity(double wavelength) const
{
    if (!_wavelengthRange.containsFuzzy(wavelength)) return 0.;

    Array params;
    double sum = 0.;
    int M = _snapshot->numEntities();
    for (int m = 0; m != M; ++m)
    {
        _snapshot->parameters(m, params);
        sum += _sedFamily->specificLuminosity(wavelength, params);
    }
    return sum;
}

////////////////////////////////////////////////////////////////////

double ImportedSource::specificLuminosity(double wavelength, int m) const
{
    if (!_wavelengthRange.containsFuzzy(wavelength)) return 0.;

    Array params;
    _snapshot->parameters(m, params);
    return _sedFamily->specificLuminosity(wavelength, params);
}

////////////////////////////////////////////////////////////////////

double ImportedSource::meanSpecificLuminosity(Range wavelengthRange, int m) const
{
    wavelengthRange.intersect(_wavelengthRange);
    if (wavelengthRange.empty()) return 0.;

    Array params, lambdav, pv, Pv;
    _snapshot->parameters(m, params);
    return _sedFamily->cdf(lambdav, pv, Pv, wavelengthRange, params) / wavelengthRange.width();
}

////////////////////////////////////////////////////////////////////

double ImportedSource::meanSpecificLuminosity(const Band* band, int m) const
{
    Range wavelengthRange = band->wavelengthRange();
    wavelengthRange.intersect(_wavelengthRange);
    if (wavelengthRange.empty()) return 0.;

    Array params, lambdav, pv, Pv;
    _snapshot->parameters(m, params);
    double Ltot = _sedFamily->cdf(lambdav, pv, Pv, wavelengthRange, params);
    return band->meanSpecificLuminosity(lambdav, pv) * Ltot;
}

////////////////////////////////////////////////////////////////////

void ImportedSource::prepareForLaunch(double sourceBias, size_t firstIndex, size_t numIndices)
{
    // skip preparation if there are no entities
    int M = _snapshot->numEntities();
    if (!M) return;

    // calculate the launch weight for each entity, normalized to unity
    _Wv = (1 - sourceBias) * _Lv + sourceBias / M;

    // determine the first history index for each entity
    _Iv.resize(M + 1);
    _Iv[0] = firstIndex;
    double W = 0.;
    for (int m = 1; m != M; ++m)
    {
        // track the cumulative normalized weight as a floating point number
        // and limit the index to firstIndex+numIndices to avoid issues with rounding errors
        W += _Wv[m - 1];
        _Iv[m] = firstIndex + min(numIndices, static_cast<size_t>(std::round(W * numIndices)));
    }
    _Iv[M] = firstIndex + numIndices;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // an instance of this class holds the normalized regular and cumulative spectral distributions for a single entity
    // which can be used to generate wavelengths and to calculate bias weights
    class EntitySED
    {
    private:
        // these two variables unambiguously identify a particular entity, even with multiple imported sources
        int _m{-1};                          // entity index
        const Snapshot* _snapshot{nullptr};  // snapshot
        Array _lambdav, _pv, _Pv;            // normalized distributions

    public:
        EntitySED() {}

        // sets the normalized distributions from the SED family if this is a different entity or snapshot
        void setIfNeeded(int m, const Snapshot* snapshot, const SEDFamily* family, Range range)
        {
            if (m != _m || snapshot != _snapshot)
            {
                Array params;
                snapshot->parameters(m, params);
                family->cdf(_lambdav, _pv, _Pv, range, params);
                _snapshot = snapshot;
                _m = m;
            }
        }

        // returns a random wavelength generated from the distribution
        double generateWavelength(Random* random) const { return random->cdfLogLog(_lambdav, _pv, _Pv); }

        // returns the normalized specific luminosity for the given wavelength
        double specificLuminosity(double lambda) const
        {
            return NR::value<NR::interpolateLogLog>(lambda, _lambdav, _pv);
        }
    };

    // setup an SED instance for each parallel execution thread to cache discretized SED data; this works even if
    // there are multiple sources of this type because each thread handles a single photon packet at a time
    thread_local EntitySED t_sed;
}

namespace
{
    // an instance of this class offers the velocity interface for an imported entity
    class EntityVelocity : public VelocityInterface
    {
    private:
        Vec _bfv;

    public:
        EntityVelocity() {}
        void setBulkVelocity(Vec bfv) { _bfv = bfv; }
        void applyVelocityDispersion(Random* random, double sigma) { _bfv += sigma * random->maxwell(); }
        Vec velocity() const override { return _bfv; }
    };

    // setup a velocity instance (with the redshift interface) for each parallel execution thread; this works even if
    // there are multiple sources of this type because each thread handles a single photon packet at a time
    thread_local EntityVelocity t_velocity;
}

////////////////////////////////////////////////////////////////////

void ImportedSource::launch(PhotonPacket* pp, size_t historyIndex, double L) const
{
    // select the entity corresponding to this history index
    auto m = std::upper_bound(_Iv.cbegin(), _Iv.cend(), historyIndex) - _Iv.cbegin() - 1;

    // if there are no entities in the source, or the selected entity has no contribution,
    // launch a photon packet with zero luminosity
    if (m < 0 || !_Lv[m])
    {
        pp->launch(historyIndex, _arbitaryWavelength, 0., Position(), Direction());
        return;
    }

    // calculate the weight related to biased source selection
    double ws = _Lv[m] / _Wv[m];

    // get the normalized regular and cumulative distributions for this entity, if not already available
    t_sed.setIfNeeded(m, _snapshot, _sedFamily, _wavelengthRange);

    // generate a random wavelength from the SED and/or from the bias distribution
    double lambda, w;
    if (!_xi)
    {
        // no biasing -- simply use the intrinsic spectral distribution
        lambda = t_sed.generateWavelength(random());
        w = 1.;
    }
    else
    {
        // biasing -- use one or the other distribution
        if (random()->uniform() > _xi)
            lambda = t_sed.generateWavelength(random());
        else
            lambda = _biasDistribution->generateWavelength();

        // calculate the compensating weight factor
        double s = t_sed.specificLuminosity(lambda);
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
            double b = _biasDistribution->probability(lambda);
            w = s / ((1 - _xi) * s + _xi * b);
        }
    }

    // generate a random position for this entity
    Position bfr = _snapshot->generatePosition(m);

    // provide a redshift interface for the appropriate velocity, if enabled
    VelocityInterface* bvi = nullptr;
    if (!_oligochromatic && _importVelocity)
    {
        t_velocity.setBulkVelocity(_snapshot->velocity(m));
        if (_importVelocityDispersion) t_velocity.applyVelocityDispersion(random(), _snapshot->velocityDispersion(m));
        bvi = &t_velocity;
    }

    // launch the photon packet with isotropic direction
    pp->launch(historyIndex, lambda, L * w * ws, bfr, random()->direction(), bvi);
}

////////////////////////////////////////////////////////////////////

const Snapshot* ImportedSource::snapshot() const
{
    return _snapshot;
}

////////////////////////////////////////////////////////////////////
