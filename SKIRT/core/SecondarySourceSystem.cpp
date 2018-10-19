/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SecondarySourceSystem.hpp"
#include "BulkVelocityInterface.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "MediumSystem.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PhotonPacket.hpp"
#include "ProbePhotonPacketInterface.hpp"
#include "ProcessManager.hpp"
#include "Random.hpp"
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

bool SecondarySourceSystem::prepareForlaunch(size_t numPackets)
{
    // calculate the absorbed (and thus to be emitted) dust luminosity for each spatial cell
    // this can be somewhat time-consuming, so we do this in parallel
    int numCells = _ms->numCells();
    _Lv.resize(numCells);
    find<ParallelFactory>()->parallelDistributed()->call(numCells, [this](size_t firstIndex, size_t numIndices)
    {
        for (size_t m=firstIndex; m!=firstIndex+numIndices; ++m)
        {
            _Lv[m] = _ms->absorbedLuminosity(m, MaterialMix::MaterialType::Dust);
        }
    });
    ProcessManager::sumToAll(_Lv);

    // calculate the total luminosity; if it is zero, report failure
    _L = _Lv.sum();
    if (_L <= 0.) return false;

    // calculate the average luminosity contribution for each packet
    _Lpp = _L / numPackets;

    // normalize the individual luminosities to unity
    _Lv /= _L;

    // determine a uniform weight for each cell, skipping zero-luminosity cells, and normalize to unity
    Array wv(numCells);
    for (int m=0; m!=numCells; ++m) wv[m] = _Lv[m] > 0 ? 1. : 0.;
    wv /= wv.sum();

    // calculate the final, composite-biased launch weight for each cell, normalized to unity
    double xi = _config->secondarySpatialBias();
    _Wv = (1-xi)*_Lv + xi*wv;

    // determine the first history index for each cell
    _Iv.resize(numCells+1);
    _Iv[0] = 0;
    for (int m=1; m!=numCells; ++m)
    {
        // limit first index to numPackets to avoid run-over due to rounding errors
        _Iv[m] = min(numPackets, _Iv[m-1] + static_cast<size_t>(std::round(_Wv[m-1] * numPackets)));
    }
    _Iv[numCells] = numPackets;

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
    class DustCellEmission : public BulkVelocityInterface
    {
    private:
        int _m{-1};                 // spatial cell index
        Array _lambdav, _pv, _Pv;   // normalized emission spectra
        Vec _bfv;                   // bulk velocity

    public:
        DustCellEmission() { }

        // calculates the emission information for the given cell if it is a different cell of what's already stored
        void calculateIfNeeded(int m, MediumSystem* ms, Configuration* config)
        {
            if (m!=_m)
            {
                _m=m;

                // copy the wavelengths of the configured emission wavelength grid
                auto wavelengthGrid = config->dustEmissionWLG();
                int numWavelengths = wavelengthGrid->numBins();
                for (int ell=0; ell!=numWavelengths; ++ell) _lambdav[ell] = wavelengthGrid->wavelength(ell);

                // get the mean intensities of the radiation field in the cell
                Array Jv = ms->meanIntensity(m);

                // accumulate the emmissivity spectrum for all dust medium components in the cell, weighed by density
                // XXXX

                // calculate the emission spectrum
                (void)config; // XXXX

                // remember the average bulk velocity
                _bfv = ms->bulkVelocity(m);
            }
        }

        // returns a random wavelength generated from the spectral distribution
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

        Vec bulkVelocity() const override { return _bfv; }
    };

    // setup a DustCellEmission instance for each parallel execution thread to cache dust emission information
    thread_local DustCellEmission t_dustcell;
}

////////////////////////////////////////////////////////////////////

void SecondarySourceSystem::launch(PhotonPacket* pp, size_t historyIndex) const
{
    // select the spatial cell from which to launch based on the history index of this photon packet
    auto m = std::upper_bound(_Iv.cbegin(), _Iv.cend(), historyIndex) - _Iv.cbegin() - 1;
    double L = _Lpp * _Lv[m] / _Wv[m];

    // calculate the emission spectrum and bulk velocity for this cell, if not already available
    t_dustcell.calculateIfNeeded(m, _ms, _config);

    // generate a random wavelength from the emission spectrum for the cell and/or from the bias distribution
    double lambda, w;
    {
        double xi = _config->secondaryWavelengthBias();
        if (!xi)
        {
            // no biasing -- simply use the intrinsic spectral distribution
            lambda = t_dustcell.generateWavelength(_random);
            w = 1.;
        }
        else
        {
            // biasing -- use one or the other distribution
            if (_random->uniform() > xi) lambda = t_dustcell.generateWavelength(_random);
            else lambda = _config->secondaryWavelengthBiasDistribution()->generateWavelength();

            // calculate the compensating weight factor
            double s = t_dustcell.specificLuminosity(lambda);
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
                double b = _config->secondaryWavelengthBiasDistribution()->probability(lambda);
                w = s / ((1-xi)*s + xi*b);
            }
        }
    }

    // generate a random position in this spatial cell
    Position bfr = _ms->grid()->randomPositionInCell(m);

    // provide a redshift interface for the appropriate velocity, if it is nonzero
    BulkVelocityInterface* bvi = t_dustcell.bulkVelocity().isNull() ? nullptr : &t_dustcell;

    // launch the photon packet with isotropic direction
    pp->launch(historyIndex, lambda, L*w, bfr, _random->direction(), bvi);

    // add origin info (we combine all medium components, so we cannot differentiate them here)
    pp->setSecondaryOrigin(0);

    // invoke launch call-back if installed
    if (_callback) _callback->probePhotonPacket(pp);
}

////////////////////////////////////////////////////////////////////
