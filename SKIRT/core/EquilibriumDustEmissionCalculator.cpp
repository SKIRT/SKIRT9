/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "EquilibriumDustEmissionCalculator.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PlanckFunction.hpp"
#include "ProcessManager.hpp"

////////////////////////////////////////////////////////////////////

void EquilibriumDustEmissionCalculator::precalculate(SimulationItem* item, const Array& lambdav, const Array& sigmaabsv)
{
    // perform initialization that needs to happen only once
    if (_rfsigmaabsvv.empty())
    {
        auto config = item->find<Configuration>();

        // obtain the simulation's radiation field wavelength grid
        auto radiationFieldWLG = config->radiationFieldWLG();
        radiationFieldWLG->setup();
        int n = radiationFieldWLG->numBins();
        _rflambdav.resize(n);
        _rfdlambdav.resize(n);
        for (int k = 0; k != n; ++k)
        {
            _rflambdav[k] = radiationFieldWLG->wavelength(k);
            _rfdlambdav[k] = radiationFieldWLG->effectiveWidth(k);
        }

        // if requested by the configuration, precalculate the CMB source term;
        // if not, the array remains initialized to all zeroes so including it has no effect
        _Bcmbv.resize(n);
        if (config->includeHeatingByCMB())
        {
            PlanckFunction B(Constants::Tcmb() * (1. + config->redshift()));
            for (int k = 0; k != n; ++k) _Bcmbv[k] = B(_rflambdav[k]);
        }

        // obtain the simulation's dust emission wavelength grid, if there is one
        if (config->hasDustEmission())
        {
            auto dustEmissionWLG = config->dustEmissionWLG();
            dustEmissionWLG->setup();
            _emlambdav = dustEmissionWLG->extlambdav();
        }

        // build the temperature grid on which we store the Planck-integrated absorption cross sections
        NR::buildPowerLawGrid(_Tv, 0., 5000., 1000, 500.);
    }

    // interpolate the absorption cross sections on the radiation field wavelength grid
    _rfsigmaabsvv.emplace_back(NR::resample<NR::interpolateLogLog>(_rflambdav, lambdav, sigmaabsv));

    // interpolate the absorption cross sections on the dust emission field wavelength grid, if there is one
    if (_emlambdav.size())
    {
        _emsigmaabsvv.emplace_back(NR::resample<NR::interpolateLogLog>(_emlambdav, lambdav, sigmaabsv));
    }

    // calculate the Planck-integrated absorption cross sections on the temperature grid
    // this can take a few seconds for all populations/size bins combined,
    // so we parallelize the loop but there is no reason to log progress
    size_t numT = _Tv.size();
    Array planckabsv(numT);
    item->find<ParallelFactory>()->parallelDistributed()->call(numT, [this, &lambdav, &sigmaabsv, &planckabsv](
                                                                         size_t firstIndex, size_t numIndices) {
        size_t numLambda = lambdav.size();
        for (size_t p = firstIndex; p != firstIndex + numIndices; ++p)
        {
            if (p)  // skip p=0, leaving the corresponding value at zero
            {
                PlanckFunction B(_Tv[p]);
                double planckabs = 0.;
                for (size_t j = 1; j != numLambda; ++j)  // skip the first wavelength so we can determine a bin width
                {
                    double lambda = lambdav[j];
                    double dlambda = lambdav[j] - lambdav[j - 1];
                    planckabs += sigmaabsv[j] * B(lambda) * dlambda;
                }
                planckabsv[p] = planckabs;
            }
        }
    });
    ProcessManager::sumToAll(planckabsv);
    _planckabsvv.emplace_back(std::move(planckabsv));
}

////////////////////////////////////////////////////////////////////

size_t EquilibriumDustEmissionCalculator::allocatedBytes() const
{
    size_t allocatedSize = 0;
    allocatedSize += _rflambdav.size();
    allocatedSize += _rfdlambdav.size();
    allocatedSize += _Bcmbv.size();
    allocatedSize += _emlambdav.size();
    allocatedSize += _Tv.size();
    if (!_rfsigmaabsvv.empty()) allocatedSize += _rfsigmaabsvv.size() * _rfsigmaabsvv[0].size();
    if (!_emsigmaabsvv.empty()) allocatedSize += _emsigmaabsvv.size() * _emsigmaabsvv[0].size();
    if (!_planckabsvv.empty()) allocatedSize += _planckabsvv.size() * _planckabsvv[0].size();
    return allocatedSize * sizeof(double);
}

////////////////////////////////////////////////////////////////////

int EquilibriumDustEmissionCalculator::numBins() const
{
    return _rfsigmaabsvv.size();
}

////////////////////////////////////////////////////////////////////

double EquilibriumDustEmissionCalculator::equilibriumTemperature(int b, const Array& Jv) const
{
    // integrate the input side of the energy balance equation
    double inputabs = (_rfsigmaabsvv[b] * (Jv + _Bcmbv) * _rfdlambdav).sum();

    // find the temperature corresponding to this amount of emission on the output side of the equation
    if (inputabs > 0.)
        return NR::clampedValue<NR::interpolateLinLin>(inputabs, _planckabsvv[b], _Tv);
    else
        return 0.;
}

////////////////////////////////////////////////////////////////////

Array EquilibriumDustEmissionCalculator::emissivity(const Array& Jv) const
{
    int numWavelengths = _emlambdav.size();
    int numBins = _rfsigmaabsvv.size();

    Array ev(numWavelengths);
    for (int b = 0; b != numBins; ++b)
    {
        double T = equilibriumTemperature(b, Jv);
        PlanckFunction B(T);
        for (int ell = 0; ell < numWavelengths; ell++)
        {
            ev[ell] += _emsigmaabsvv[b][ell] * B(_emlambdav[ell]);
        }
    }
    return ev;
}

////////////////////////////////////////////////////////////////////
