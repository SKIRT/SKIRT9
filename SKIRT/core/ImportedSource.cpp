/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedSource.hpp"
#include "Constants.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PhotonPacket.hpp"
#include "ProcessManager.hpp"
#include "Random.hpp"
#include "RedshiftInterface.hpp"
#include "SEDFamily.hpp"
#include "Snapshot.hpp"
#include "WavelengthRangeInterface.hpp"

////////////////////////////////////////////////////////////////////

void ImportedSource::setupSelfAfter()
{
    Source::setupSelfAfter();

    // get the primary source wavelength range
    _wavelengthRange = interface<WavelengthRangeInterface>()->wavelengthRange();

    // create the snapshot with preconfigured spatial columns
    _snapshot = createAndOpenSnapshot();

    // add optional columns if applicable
    if (_importVelocity) _snapshot->importVelocity();
    _snapshot->importParameters(_sedFamily->parameterInfo());

    // read the data from file
    _snapshot->readAndClose();

    // construct a vector with the (normalized) luminosity for each entity
    int M = _snapshot->numEntities();
    if (M)
    {
        find<Log>()->info("Calculating luminosities from imported properties");
        _Lv.resize(M);

        // integrating over the SED for each entity can be time-consuming, so we do this in parallel
        find<ParallelFactory>()->parallelDistributed()->call(M, [this](size_t firstIndex, size_t numIndices)
        {
            Array lambdav, pv, Pv;  // the contents of these arrays is not used, so this could be optimized if needed
            Array params;

            for (size_t m=firstIndex; m!=firstIndex+numIndices; ++m)
            {
                _snapshot->parameters(m, params);
                _Lv[m] = _sedFamily->cdf(lambdav, pv, Pv, _wavelengthRange, params);
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
    for (int m=0; m!=M; ++m)
    {
        _snapshot->parameters(m, params);
        sum += _sedFamily->specificLuminosity(wavelength, params);
    }
    return sum;
}

////////////////////////////////////////////////////////////////////

void ImportedSource::prepareForLaunch(double sourceBias, size_t firstIndex, size_t numIndices)
{
    // skip preparation if there are no entities
    int M = _snapshot->numEntities();
    if (!M) return;

    // calculate the launch weight for each entity, normalized to unity
    _Wv = (1-sourceBias)*_Lv + sourceBias/M;

    // determine the first history index for each entity
    _Iv.resize(M+1);
    _Iv[0] = firstIndex;
    for (int m=1; m!=M; ++m)
    {
        // limit first index to last index to avoid run-over due to rounding errors
        _Iv[m] = min(firstIndex+numIndices, _Iv[m-1] + static_cast<size_t>(std::round(_Wv[m-1] * numIndices)));
    }
    _Iv[M] = firstIndex+numIndices;
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
        int _m{-1};                         // entity index
        const Snapshot* _snapshot{nullptr}; // snapshot
        Array _lambdav, _pv, _Pv;           // normalized distributions

    public:
        EntitySED() { }

        // sets the normalized distributions from the SED family if this is a different entity or snapshot
        void setIfNeeded(int m, const Snapshot* snapshot, const SEDFamily* family, Range range)
        {
            if (m!=_m || snapshot!=_snapshot)
            {
                Array params;
                snapshot->parameters(m, params);
                family->cdf(_lambdav, _pv, _Pv, range, params);
                _snapshot=snapshot; _m=m;
            }
        }

        // returns a random wavelength generated from the distribution
        double generateWavelength(Random* random) const
        {
            return random->cdfLogLog(_lambdav, _pv, _Pv);
        }

        // returns the normalized specific luminosity for the given wavelength
        double specificLuminosity(double lambda) const
        {
            int i = NR::locateClip(_lambdav, lambda);
            return NR::interpolateLogLog(lambda, _lambdav[i], _lambdav[i+1], _pv[i], _pv[i+1]);
        }
    };

    // setup an SED instance for each parallel execution thread to cache discretized SED data; this works even if
    // there are multiple sources of this type because each thread handles a single photon packet at a time
    thread_local EntitySED t_sed;
}

namespace
{
    // an instance of this class offers the redshift interface for the velocity specified for an imported entity
    class EntityVelocity : public RedshiftInterface
    {
    private:
        Vec _bfv;
    public:
        EntityVelocity() { }
        void setVelocity(Vec bfv) { _bfv = bfv; }
        double redshiftForDirection(Direction bfk) const override { return -Vec::dot(_bfv, bfk) / Constants::c(); }
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
        pp->launch(historyIndex, _wavelengthRange.mid(), 0., Position(), Direction());
        return;
    }

    // calculate the weight related to biased source selection
    double ws = _Lv[m] / _Wv[m];

    // get the normalized regular and cumulative distributions for this entity, if not already available
    t_sed.setIfNeeded(m, _snapshot, _sedFamily, _wavelengthRange);

    // generate a random wavelength from the SED and/or from the bias distribution
    double lambda, w;
    double xi = wavelengthBias();
    if (!xi)
    {
        // no biasing -- simply use the intrinsic spectral distribution
        lambda = t_sed.generateWavelength(random());
        w = 1.;
    }
    else
    {
        // biasing -- use one or the other distribution
        if (random()->uniform() > xi) lambda = t_sed.generateWavelength(random());
        else lambda = wavelengthBiasDistribution()->generateWavelength();

        // calculate the compensating weight factor
        double s = t_sed.specificLuminosity(lambda);
        if (!s)
        {
            // if the wavelength can't occur in the intrinsic distribution,
            // the weight factor is zero regardless of the probability in the bias distribution
            // (handling this separately also avoids NaNs in the pathological case s=d=0)
            w = 0.;
        }
        else
        {
            // regular composite bias weight
            double b = wavelengthBiasDistribution()->probability(lambda);
            w = s / ((1-xi)*s + xi*b);
        }
    }

    // generate a random position for this entity
    Position bfr = _snapshot->generatePosition(m);

    // provide a redshift interface for the appropriate velocity, if enabled
    RedshiftInterface* rsi = nullptr;
    if (importVelocity())
    {
        t_velocity.setVelocity(_snapshot->velocity(m));
        rsi = &t_velocity;
    }

    // launch the photon packet with isotropic direction
    pp->launch(historyIndex, lambda, L*w*ws, bfr, random()->direction(), rsi);
}

////////////////////////////////////////////////////////////////////
