/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LineGasSecondarySource.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "EmittingGasMix.hpp"
#include "FatalError.hpp"
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

////////////////////////////////////////////////////////////////////

LineGasSecondarySource::LineGasSecondarySource(SimulationItem* parent, int h) : SecondarySource(parent), _h(h) {}

////////////////////////////////////////////////////////////////////

double LineGasSecondarySource::prepareLuminosities()
{
    // cache some pointers for later use
    _config = find<Configuration>();
    _ms = find<MediumSystem>();
    _mix = dynamic_cast<const EmittingGasMix*>(_ms->mix(0, _h));
    _random = find<Random>();

    // copy the line wavelengths into local storage
    _centers = _mix->lineEmissionCenters();
    _numLines = _centers.size();
    if (!_numLines) throw FATALERROR("Line emitting gas mix does not declare any lines");

    // determine whether the component's medium state has a temperature variable
    for (const StateVariable& variable : _mix->specificStateVariableInfo())
        if (variable.identifier() == StateVariable::Identifier::Temperature) _hasTemperature = true;

    // if so, also retrieve the particle masses corresponding to each line
    if (_hasTemperature) _masses = _mix->lineEmissionMasses();

    // calculate the bolometric luminosity for each spatial cell by adding all line luminosities
    // this can be time-consuming, so we do this in parallel
    int numCells = _ms->numCells();
    _Lv.resize(numCells);
    find<ParallelFactory>()->parallelDistributed()->call(numCells, [this](size_t firstIndex, size_t numIndices) {
        for (size_t m = firstIndex; m != firstIndex + numIndices; ++m)
        {
            _Lv[m] = _ms->lineEmissionSpectrum(m, _h).sum();
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
    log->info("Line luminosity for gas medium component " + std::to_string(_h) + ": "
              + StringUtils::toString(units->obolluminosity(L), 'g') + " " + units->ubolluminosity());
    log->info("  Emitting from " + std::to_string(emittingCells) + " out of " + std::to_string(numCells)
              + " spatial cells");

    // return the total luminosity
    return L;
}

////////////////////////////////////////////////////////////////////

void LineGasSecondarySource::preparePacketMap(size_t firstIndex, size_t numIndices)
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
    // from the line spectrum of the gas component being handled in a given cell in the spatial grid,
    // and remembers the information for fast retrieval. This information includes the normalized regular
    // and cumulative discrete distributions of the line luminosities, calculated from the stored radiation
    // field and the gas material properties, and the average bulk velocity of the material in the cell,
    // obtained from the medium system.
    class GasCellEmission : public VelocityInterface
    {
    private:
        // information on a particular spatial cell, initialized by calculateIfNeeded()
        int _m{-1};      // spatial cell index
        Array _pv, _Pv;  // normalized discrete emission spectrum (relative line luminosities)
        Vec _vbulk;      // bulk velocity

        // information on a particular photon packet launch, initialized by generateThermalVelocity()
        Vec _vtherm;  // random thermal velocity

    public:
        // instances of this class are allocated in thread-local storage, which means that the constructor is
        // invoked under the hood for each new execution thread; hence it does trivial initialization only
        GasCellEmission() {}

        // calculates the emission information for the given cell if it has not already been calculated
        //   m:  cell index
        //   h:  medium component index
        //   ms: medium system
        void calculateIfNeeded(int m, int h, MediumSystem* ms)
        {
            // if this photon packet is launched from the same cell as the previous one, we don't need to do anything
            if (m == _m) return;

            // remember the new cell index
            _m = m;

            // get the emission spectrum of this medium component in the cell
            _pv = ms->lineEmissionSpectrum(m, h);

            // calculate the normalized plain and cumulative distributions
            double norm = NR::cdf(_Pv, _pv);
            _pv /= norm;

            // get the average bulk velocity for this cell
            _vbulk = ms->bulkVelocity(m);

            // reset the thermal velocity
            _vtherm = Vec();
        }

        // returns a random line index generated from the discrete spectral distribution
        int generateLineIndex(Random* random) const { return NR::locateClip(_Pv, random->uniform()); }

        // returns the normalized line luminosity for the given wavelength
        double luminosity(int index) const { return _pv[index]; }

        // generates and stores a random thermal velocity for the given medium temperature and particle mass
        void generateThermalVelocity(Random* random, double T, double M)
        {
            if (T > 0. && M > 0.)
                _vtherm = sqrt(Constants::k() * T / M) * random->maxwell();
            else
                _vtherm = Vec();
        }

        // returns the sum of the bulk velocity and the random thermal velocity
        Vec velocity() const override { return _vbulk + _vtherm; }
    };

    // setup an instance of the above class to cache emission information for each parallel execution thread
    thread_local GasCellEmission t_gascell;
}

////////////////////////////////////////////////////////////////////

void LineGasSecondarySource::launch(PhotonPacket* pp, size_t historyIndex, double L) const
{
    // select the spatial cell from which to launch based on the history index of this photon packet
    auto m = std::upper_bound(_Iv.cbegin(), _Iv.cend(), historyIndex) - _Iv.cbegin() - 1;

    // calculate the weight related to biased source selection
    double ws = _Lv[m] / _Wv[m];

    // calculate the emission spectrum and bulk velocity for this cell, if not already available
    t_gascell.calculateIfNeeded(m, _h, _ms);

    // randomly select a line from the luminosity distribution and/or from a uniform distribution
    int index = 0;
    double w = 1.;
    if (_numLines > 1)
    {
        double xi = _mix->wavelengthBias();
        if (!xi)
        {
            // no biasing -- simply use the intrinsic spectral distribution
            index = t_gascell.generateLineIndex(_random);
        }
        else
        {
            // biasing -- use one or the other distribution
            if (_random->uniform() > xi)
                index = t_gascell.generateLineIndex(_random);
            else
                index = _random->uniform() * _numLines;

            // calculate the compensating composite bias weight factor
            double s = t_gascell.luminosity(index);
            double b = 1. / _numLines;
            w = s / ((1 - xi) * s + xi * b);
        }
    }

    // if a temperature is available, generate a random thermal velocity
    if (_hasTemperature) t_gascell.generateThermalVelocity(_random, _ms->temperature(m, _h), _masses[index]);

    // generate a random position in this spatial cell
    Position bfr = _ms->grid()->randomPositionInCell(m);

    // provide a redshift interface for the appropriate velocity, if it is nonzero
    VelocityInterface* bvi = t_gascell.velocity().isNull() ? nullptr : &t_gascell;

    // launch the photon packet with isotropic direction
    pp->launch(historyIndex, _centers[index], L * ws * w, bfr, _random->direction(), bvi);
}

////////////////////////////////////////////////////////////////////
