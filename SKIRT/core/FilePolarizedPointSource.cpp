/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FilePolarizedPointSource.hpp"
#include "AngularDistributionInterface.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PhotonPacket.hpp"
#include "PolarizationProfileInterface.hpp"
#include "ProcessManager.hpp"
#include "Random.hpp"
#include "TabulatedSED.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // this helper class integrates the intensity tabulated in the user input over the unit sphere
    // and exposes the resulting mean spectrum (with arbitrary normalization) as a regular ContSED subclass.
    class MeanSED : public TabulatedSED
    {
    public:
        // this constructor remembers a reference to the intensity table, installs the newly created object
        // as a child to the specified parent so it will automatically be deleted, and calls its setup() function
        explicit MeanSED(SimulationItem* parent, const StoredTable<2>& tableI) : _tableI(tableI)
        {
            parent->addChild(this);
            setup();
        }

    protected:
        // this function calculates the unit-sphere-integrated spectrum and stores it into the provided argument arrays
        void getWavelengthsAndLuminosities(Array& lambdav, Array& pv) const override
        {
            // get the cosine grid points and the distances between the grid points
            Array cosv;
            _tableI.axisArray<0>(cosv);
            size_t numCos = cosv.size();
            Array deltacosv(numCos);
            for (size_t t = 1; t != numCos; ++t) deltacosv = cosv[t] - cosv[t - 1];

            // get the wavelength grid points
            _tableI.axisArray<1>(lambdav);
            size_t numLambda = lambdav.size();

            // loop over all wavelengths, in parallel
            pv.resize(numLambda);
            find<ParallelFactory>()->parallelDistributed()->call(numLambda, [this, numCos, &deltacosv, &pv](
                                                                                size_t firstIndex, size_t numIndices) {
                for (size_t ell = firstIndex; ell != firstIndex + numIndices; ++ell)
                {
                    // perform the integration using the trapezium rule, with arbitrary normalization
                    double total = 0.;
                    for (size_t t = 1; t != numCos; ++t)
                    {
                        total += (_tableI.valueAtIndices(t - 1, ell) + _tableI.valueAtIndices(t, ell)) * deltacosv[t];
                    }
                    pv[ell] = total;
                }
            });
            ProcessManager::sumToAll(pv);
        }

    private:
        // a reference to the intensity table
        const StoredTable<2>& _tableI;
    };
}

////////////////////////////////////////////////////////////////////

void FilePolarizedPointSource::setupSelfBefore()
{
    Source::setupSelfBefore();

    // open the tables with Stokes vector components as a function of wavelength and inclination angle cosine
    _tableI.open(this, filename(), "costheta(1),lambda(m)", "I(W/m)", false, false);
    _tableQ.open(this, filename(), "costheta(1),lambda(m)", "Q(W/m)", false, false);
    _tableU.open(this, filename(), "costheta(1),lambda(m)", "U(W/m)", false, false);
    _tableV.open(this, filename(), "costheta(1),lambda(m)", "V(W/m)", false, false);

    // construct the mean SED (averaged over the unit sphere) from the input file
    _sed = new MeanSED(this, _tableI);

    // get the symmetry axis orientation
    Vec sym(symmetryX(), symmetryY(), symmetryZ());
    double norm = sym.norm();
    if (norm == 0) throw FATALERROR("Symmetry axis direction cannot be null vector");
    _sym = Direction(sym / norm);

    // cache wavelength sampling parameters depending on whether this is an oligo- or panchromatic simulation
    auto config = find<Configuration>();
    if (config->oligochromatic())
    {
        _xi = 1.;  // always use bias distribution for oligochromatic simulations
        _biasDistribution = config->oligoWavelengthBiasDistribution();
    }
    else
    {
        _xi = wavelengthBias();
        _biasDistribution = wavelengthBiasDistribution();
    }

    // if we have a nonzero bulk velocity, set the interface object to ourselves; otherwise leave it at null pointer
    if (hasVelocity()) _bvi = this;
}

////////////////////////////////////////////////////////////////////

void FilePolarizedPointSource::setupSelfAfter()
{
    Source::setupSelfAfter();

    // warn the user if this source's intrinsic wavelength range does not fully cover the configured wavelength range
    informAvailableWavelengthRange(_sed->intrinsicWavelengthRange(), _sed->type());
}

////////////////////////////////////////////////////////////////////

int FilePolarizedPointSource::dimension() const
{
    if (positionX() || positionY() || symmetryX() || symmetryY()) return 3;
    if (hasVelocity() && (velocityX() || velocityY())) return 3;
    return 2;
}

////////////////////////////////////////////////////////////////////

bool FilePolarizedPointSource::hasVelocity() const
{
    if (velocityX() || velocityY() || velocityZ())
    {
        // refuse velocity for oligochromatic simulations
        // (this function is called from Configure so we cannot precompute this during setup)
        auto config = find<Configuration>();
        if (!config->oligochromatic()) return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////

Vec FilePolarizedPointSource::velocity() const
{
    return Vec(velocityX(), velocityY(), velocityZ());
}

////////////////////////////////////////////////////////////////////

Range FilePolarizedPointSource::wavelengthRange() const
{
    return _sed->normalizationWavelengthRange();
}

////////////////////////////////////////////////////////////////////

double FilePolarizedPointSource::luminosity() const
{
    return _normalization->luminosityForSED(_sed);
}

////////////////////////////////////////////////////////////////////

double FilePolarizedPointSource::specificLuminosity(double wavelength) const
{
    return _sed->normalizationWavelengthRange().containsFuzzy(wavelength)
               ? _sed->specificLuminosity(wavelength) * luminosity()
               : 0.;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // this helper class represents the angular distribution at a given wavelength,
    // derived from the intensity tabulated in the user input as a function of inclination,
    // and taking into account the orientation of the symmetry axis
    class LocalAngularDistribution : public AngularDistributionInterface
    {
    public:
        // this function obtains the angular probability distribution for the specified wavelength
        // and stores the symmetry axis orientation and the simulation's random generator
        void initialize(const StoredTable<2>& tableI, double lambda, Direction sym, Random* random)
        {
            tableI.cdf(_cosv, _pv, _Pv, Range(-1., 1.), lambda);
            _sym = sym;
            _random = random;
        }

        // this function returns the probability for the angle between the given direction and the symmetry axis
        double probabilityForDirection(Direction bfk) const override
        {
            double cosine = Vec::dot(_sym, bfk);
            return NR::value<NR::interpolateLinLin>(cosine, _cosv, _pv);
        }

        // this function generates a random cosine from the probability distribution, and then
        // returns a direction with the corresponding angle relative to the symmetry axis and a random azimuth
        Direction generateDirection() const
        {
            double cosine = _random->cdfLinLin(_cosv, _Pv);
            return _random->direction(_sym, cosine);
        }

    private:
        Array _cosv;  // cosine grid points
        Array _pv;    // corresponding normalized probability distribution
        Array _Pv;    // corresponding normalized cumulative probability distribution

        Direction _sym;   // orientation of the symmetry axis
        Random* _random;  // the simulation's random generator
    };

    // setup an angular distribution instance for each parallel execution thread; this works even if
    // there are multiple sources of this type because each thread handles a single photon packet at a time
    thread_local LocalAngularDistribution t_angularDistribution;
}

////////////////////////////////////////////////////////////////////

void FilePolarizedPointSource::launch(PhotonPacket* pp, size_t historyIndex, double L) const
{
    // generate a random wavelength from the SED and/or from the bias distribution
    double lambda, w;
    if (!_xi)
    {
        // no biasing -- simply use the intrinsic spectral distribution
        lambda = _sed->generateWavelength();
        w = 1.;
    }
    else
    {
        // biasing -- use one or the other distribution
        if (random()->uniform() > _xi)
            lambda = _sed->generateWavelength();
        else
            lambda = _biasDistribution->generateWavelength();

        // calculate the compensating weight factor
        double s = _sed->specificLuminosity(lambda);
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

    // get the source position
    Position bfr(positionX(), positionY(), positionZ());

    // initialize thread-local angular distribution object for this wavelength point
    t_angularDistribution.initialize(_tableI, lambda, _sym, random());

    // get or construct polarization profile object for this wavelength point
    PolarizationProfileInterface* polarizationProfile = nullptr;

    // launch the photon packet with the proper wavelength, weight, position, direction,
    // bulk veloclity, angular distribution and polarization profile
    pp->launch(historyIndex, lambda, L * w, bfr, t_angularDistribution.generateDirection(), _bvi,
               &t_angularDistribution, polarizationProfile);
}

////////////////////////////////////////////////////////////////////
