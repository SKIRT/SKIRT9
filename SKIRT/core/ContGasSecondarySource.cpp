/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ContGasSecondarySource.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "EmittingGasMix.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PhotonPacket.hpp"
#include "ProcessManager.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "Units.hpp"
#include "VelocityInterface.hpp"
#include "WavelengthDistribution.hpp"

////////////////////////////////////////////////////////////////////

ContGasSecondarySource::ContGasSecondarySource(SimulationItem* parent, int h) : SecondarySource(parent), _h(h) {}

////////////////////////////////////////////////////////////////////

double ContGasSecondarySource::prepareLuminosities()
{
    // cache some pointers for later use
    _config = find<Configuration>();
    _ms = find<MediumSystem>();
    _mix = dynamic_cast<const EmittingGasMix*>(_ms->mix(0, _h));
    _random = find<Random>();

    // calculate the bolometric luminosity for each spatial cell by integrating the spectrum over the wavelength grid
    // this can be time-consuming, so we do this in parallel
    int numCells = _ms->numCells();
    _Lv.resize(numCells);
    find<ParallelFactory>()->parallelDistributed()->call(numCells, [this](size_t firstIndex, size_t numIndices) {
        const Array& dlambdav = _mix->emissionWavelengthGrid()->extdlambdav();
        for (size_t m = firstIndex; m != firstIndex + numIndices; ++m)
        {
            _Lv[m] = (_ms->continuumEmissionSpectrum(m, _h) * dlambdav).sum();
        }
    });
    ProcessManager::sumToAll(_Lv);

    // calculate the total luminosity, and normalize the individual luminosities to unity
    double L = _Lv.sum();
    _Lv /= L;

    // log the luminosity and the number of emitting cells
    auto log = find<Log>();
    auto units = find<Units>();
    int emittingCells = 0;  // number of nonzero luminosity cells
    for (int m = 0; m != numCells; ++m)
        if (_Lv[m] > 0.) emittingCells++;
    log->info("Continuum luminosity for gas medium component " + std::to_string(_h) + ": "
              + StringUtils::toString(units->obolluminosity(L), 'g') + " " + units->ubolluminosity());
    log->info("  Emitting from " + std::to_string(emittingCells) + " out of " + std::to_string(numCells)
              + " spatial cells");

    // return the total luminosity
    return L;
}

////////////////////////////////////////////////////////////////////

void ContGasSecondarySource::preparePacketMap(size_t firstIndex, size_t numIndices)
{
    int numCells = _ms->numCells();

    // determine a uniform weight for each cell with non-negligable emission, and normalize to unity
    Array wv(numCells);
    for (int m = 0; m != numCells; ++m) wv[m] = _Lv[m] > 0. ? 1. : 0.;
    wv /= wv.sum();

    // calculate the final, composite-biased launch weight for each cell, normalized to unity
    double xi = _config->secondarySpatialBias();
    _Wv = (1 - xi) * _Lv + xi * wv;

    // determine the first history index for each cell
    _Iv.resize(numCells + 1);
    _Iv[0] = firstIndex;
    double W = 0.;
    for (int m = 1; m != numCells; ++m)
    {
        // track the cumulative normalized weight as a floating point number
        // and limit the index to firstIndex+numIndices to avoid issues with rounding errors
        W += _Wv[m - 1];
        _Iv[m] = firstIndex + min(numIndices, static_cast<size_t>(std::round(W * numIndices)));
    }
    _Iv[numCells] = firstIndex + numIndices;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // An instance of this class obtains and/or calculates the information needed to launch photon packets
    // from the continuum spectrum of the gas component being handled in a given cell in the spatial grid,
    // and remembers the information for fast retrieval. This information includes the normalized regular
    // and cumulative emission spectra, calculated from the stored radiation field and the gas material
    // properties, and the average bulk velocity of the material in the cell, obtained from the medium system.
    class GasCellEmission : public VelocityInterface
    {
    private:
        // information initialized once, during the first call to calculateIfNeeded()
        const EmittingGasMix* _mix{nullptr};  // the gas mix;
        Array _wavelengthGrid;                // the emission wavelength grid
        Range _wavelengthRange;               // the range of the emission wavelength grid
        int _numWavelengths{0};               // the number of wavelengths in the emission wavelength grid

        // information on a particular spatial cell, initialized by calculateIfNeeded()
        int _m{-1};                // spatial cell index
        Array _lambdav, _pv, _Pv;  // normalized emission spectrum
        Vec _bfv;                  // bulk velocity

    public:
        // instances of this class are allocated in thread-local storage, which means that the constructor is
        // invoked under the hood for each new execution thread; hence it does trivial initialization only
        GasCellEmission() {}

        // calculates the emission information for the given cell if it has not already been calculated
        //   m:  cell index
        //   h:  medium component index
        //   mix: material mix for the medium component
        //   ms: medium system
        void calculateIfNeeded(int m, int h, const EmittingGasMix* mix, MediumSystem* ms)
        {
            // when called for the first time for a given gas mix, cache info on the emission wavelength grid
            if (_mix != mix)
            {
                _m = -1;
                _mix = mix;
                auto wavelengthGrid = mix->emissionWavelengthGrid();
                _wavelengthGrid = wavelengthGrid->extlambdav();
                _wavelengthRange = wavelengthGrid->wavelengthRange();
                _numWavelengths = _wavelengthGrid.size();
            }

            // if this photon packet is launched from the same cell as the previous one, we don't need to do anything
            if (m == _m) return;

            // remember the new cell index
            _m = m;

            // get the emission spectrum of this medium component in the cell
            const Array& Lv = ms->continuumEmissionSpectrum(m, h);

            // calculate the normalized plain and cumulative distributions
            NR::cdf<NR::interpolateLogLog>(_lambdav, _pv, _Pv, _wavelengthGrid, Lv, _wavelengthRange);

            // get the average bulk velocity for this cell
            _bfv = ms->bulkVelocity(m);
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

    // setup an instance of the above class to cache emission information for each parallel execution thread
    thread_local GasCellEmission t_gascell;
}

////////////////////////////////////////////////////////////////////

void ContGasSecondarySource::launch(PhotonPacket* pp, size_t historyIndex, double L) const
{
    // select the spatial cell from which to launch based on the history index of this photon packet
    auto m = std::upper_bound(_Iv.cbegin(), _Iv.cend(), historyIndex) - _Iv.cbegin() - 1;

    // calculate the weight related to biased source selection
    double ws = _Lv[m] / _Wv[m];

    // calculate the emission spectrum and bulk velocity for this cell, if not already available
    t_gascell.calculateIfNeeded(m, _h, _mix, _ms);

    // generate a random wavelength from the emission spectrum for the cell and/or from the bias distribution
    double lambda, w;
    {
        double xi = _mix->wavelengthBias();
        if (!xi)
        {
            // no biasing -- simply use the intrinsic spectral distribution
            lambda = t_gascell.generateWavelength(_random);
            w = 1.;
        }
        else
        {
            // biasing -- use one or the other distribution
            if (_random->uniform() > xi)
                lambda = t_gascell.generateWavelength(_random);
            else
                lambda = _mix->wavelengthBiasDistribution()->generateWavelength();

            // calculate the compensating weight factor
            double s = t_gascell.specificLuminosity(lambda);
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
                double b = _mix->wavelengthBiasDistribution()->probability(lambda);
                w = s / ((1 - xi) * s + xi * b);
            }
        }
    }

    // generate a random position in this spatial cell
    Position bfr = _ms->grid()->randomPositionInCell(m);

    // provide a redshift interface for the appropriate velocity, if it is nonzero
    VelocityInterface* bvi = t_gascell.velocity().isNull() ? nullptr : &t_gascell;

    // launch the photon packet with isotropic direction
    pp->launch(historyIndex, lambda, L * ws * w, bfr, _random->direction(), bvi);
}

////////////////////////////////////////////////////////////////////
