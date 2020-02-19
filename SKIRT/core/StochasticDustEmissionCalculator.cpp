/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StochasticDustEmissionCalculator.hpp"
#include "ArrayTable.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PlanckFunction.hpp"
#include "ProcessManager.hpp"
#include <unordered_map>

////////////////////////////////////////////////////////////////////

// container classes that are highly specialized to optimize the operations in this class
namespace
{
    // square matrix
    template<typename T> class Square
    {
    private:
        size_t _n;
        T* _v;

    public:
        // constructor sets maximum size and logical size
        Square(size_t n) : _n(n), _v(new T[n * n]) {}
        ~Square() { delete[] _v; }

        // sets logical size, which must not be larger than maximum size set in constructor (is not checked)
        // does not clear values; does not resize underlying memory
        void resize(size_t n) { _n = n; }

        // access to values  (const version currently not needed)
        T& operator()(size_t i, size_t j) { return _v[i * _n + j]; }
    };

    // square matrix with only items below the diagonal (i>j)
    template<typename T> class Triangle
    {
    private:
        T* _v;
        static size_t offset(size_t i) { return ((i - 1) * i) >> 1; }

    public:
        // constructor sets size (can't be changed)
        Triangle(size_t n) : _v(new T[offset(n)]) {}
        ~Triangle() { delete[] _v; }

        // access to values; must have i>j (is not checked)
        const T& operator()(size_t i, size_t j) const { return _v[offset(i) + j]; }
        T& operator()(size_t i, size_t j) { return _v[offset(i) + j]; }
    };
}

////////////////////////////////////////////////////////////////////

// helper class to construct and store a particular temperature grid and the black-body emissivity spectrum
// discretized on this temperature grid and on the output wavelength grid;
// all members are public for ease of use in SDE_Calculator
class SDE_TemperatureGrid
{
public:
    // temperature grid (indexed on p)
    Array _Tv;

    // black-body radiation on this temperature grid and on the output wavelength grid (indexed on p and ell)
    ArrayTable<2> _Bvv;

    // constructor
    SDE_TemperatureGrid(const Array& emlambdav, double Tmin, double Tmax, int NT, double ratio)
    {
        // build temperature grid (linear if ratio==1)
        NR::buildPowerLawGrid(_Tv, Tmin, Tmax, NT - 1, ratio);

        // pre-calculate black-body radiation
        int numLambda = emlambdav.size();
        _Bvv.resize(NT, numLambda);
        for (int p = 0; p != NT; ++p)
        {
            PlanckFunction B(_Tv[p]);
            for (int ell = 0; ell != numLambda; ++ell) _Bvv(p, ell) = B(emlambdav[ell]);
        }
    }

    // return estimate of allocated memory size
    size_t allocatedBytes() const { return (_Tv.size() + _Bvv.size()) * sizeof(double); }
};

////////////////////////////////////////////////////////////////////

// helper class to calculate the temperature probability distribution on a particular fixed temperature grid
// for a representative grain embedded in a given radiation field (discretized on the input wavelength grid)
// with support for pre-calculating and storing relevant data relating to the grain/temperature grid combination
class SDE_Calculator
{
private:
    const SDE_TemperatureGrid* _grid;  // temperature grid with corresponding black body discretization
    Triangle<double> _HRm;             // heating rates (indexed on f,i)
    Triangle<short> _Km;               // radiation field wavelength index k (indexed on f,i)
    Array _CRv;                        // cooling rates (indexed on p)
    Array _planckabsv;                 // Planck-integrated absorption cross sections (indexed on p)

    const Array& _rfdlambdav;  // reference to radiation field wavelength grid bin widths (indexed on k)
    const Array& _emlambdav;   // reference to dust emission wavelength grid (indexed on ell)
    Array _rfsigmaabsv;        // cross sections on radiation field wavelength grid (indexed on k)
    Array _emsigmaabsv;        // cross sections on dust emission wavelength grid (indexed on ell)

public:
    SDE_Calculator(SimulationItem* item, const SDE_TemperatureGrid* grid, const WavelengthGrid* rfWLG,
                   const Array& rflambdav, const Array& rfdlambdav, const Array& emlambdav, const Array& lambdav,
                   const Array& sigmaabsv, double bulkDensity, double meanMass, const StoredTable<1>& enthalpy)
        : _grid(grid), _HRm(grid->_Tv.size()), _Km(grid->_Tv.size()), _CRv(grid->_Tv.size()),
          _planckabsv(grid->_Tv.size()), _rfdlambdav(rfdlambdav), _emlambdav(emlambdav)
    {
        // get shortcuts to the temperature grid
        const Array& Tv = grid->_Tv;
        int NT = Tv.size();

        // calculate the enthalpy of a single dust grain in this population across the temperature grid
        Array Hv(NT);
        for (int p = 0; p < NT; p++) Hv[p] = enthalpy(Tv[p]) * meanMass / bulkDensity;

        // calculate the enthalpy bin widths
        Array dHv(NT);
        dHv[0] = Hv[1] - Hv[0];
        for (int p = 1; p < NT - 1; p++)
        {
            double Tmin = (Tv[p - 1] + Tv[p]) / 2.;
            double Tmax = (Tv[p + 1] + Tv[p]) / 2.;
            dHv[p] = (enthalpy(Tmax) - enthalpy(Tmin)) * meanMass / bulkDensity;
        }
        dHv[NT - 1] = Hv[NT - 1] - Hv[NT - 2];

        // calculate the heating rates, barring the dependency on the radiation field
        double hc = Constants::h() * Constants::c();
        for (int f = 1; f < NT; f++)
        {
            for (int i = 0; i < f; i++)
            {
                double Hdiff = Hv[f] - Hv[i];
                double lambda = hc / Hdiff;
                int k = rfWLG->bin(lambda);
                _Km(f, i) = k;
                if (k >= 0)
                {
                    double sigmaabs = NR::value<NR::interpolateLogLog>(lambda, lambdav, sigmaabsv);
                    _HRm(f, i) = hc * sigmaabs * dHv[f] / (Hdiff * Hdiff * Hdiff);
                }
            }
        }

        // calculate the cooling rates
        // this can take a few seconds for all populations/size bins combined,
        // so we parallelize the loop but there is no reason to log progress
        item->find<ParallelFactory>()->parallelDistributed()->call(
            NT, [this, &lambdav, &sigmaabsv, &Tv, &Hv](size_t firstIndex, size_t numIndices) {
                size_t numLambda = lambdav.size();
                for (size_t p = firstIndex; p != firstIndex + numIndices; ++p)
                {
                    if (p)  // skip p=0, leaving the corresponding value at zero
                    {
                        PlanckFunction B(Tv[p]);
                        double planckabs = 0.;
                        for (size_t j = 1; j != numLambda;
                             ++j)  // skip the first wavelength so we can determine a bin width
                        {
                            double lambda = lambdav[j];
                            double dlambda = lambdav[j] - lambdav[j - 1];
                            planckabs += sigmaabsv[j] * B(lambda) * dlambda;
                        }
                        _planckabsv[p] = planckabs;
                        double Hdiff = Hv[p] - Hv[p - 1];
                        _CRv[p] = planckabs / Hdiff;
                    }
                }
            });
        ProcessManager::sumToAll(_planckabsv);
        ProcessManager::sumToAll(_CRv);

        // obtain the cross sections on the input and output wavelength grids
        _rfsigmaabsv = NR::resample<NR::interpolateLogLog>(rflambdav, lambdav, sigmaabsv);
        _emsigmaabsv = NR::resample<NR::interpolateLogLog>(emlambdav, lambdav, sigmaabsv);
    }

    // return estimate of allocated memory size
    size_t allocatedBytes() const
    {
        // triangles
        size_t n = _grid->_Tv.size();
        size_t allocatedBytes = (((n - 1) * n) >> 1) * (sizeof(_HRm(0, 0)) + sizeof(_Km(0, 0)));

        // arrays
        allocatedBytes += _CRv.size() * sizeof(_CRv[0]);
        allocatedBytes += _planckabsv.size() * sizeof(_planckabsv[0]);
        allocatedBytes += _rfsigmaabsv.size() * sizeof(_rfsigmaabsv[0]);
        allocatedBytes += _emsigmaabsv.size() * sizeof(_emsigmaabsv[0]);
        return allocatedBytes;
    }

    // return the equilibrium temperature of the population
    double equilibriumTemperature(const Array& Jv) const
    {
        // integrate the input side of the energy balance equation
        double inputabs = (_rfsigmaabsv * Jv * _rfdlambdav).sum();

        // find the temperature corresponding to this amount of emission on the output side of the equation
        if (inputabs > 0.)
            return NR::clampedValue<NR::interpolateLinLin>(inputabs, _planckabsv, _grid->_Tv);
        else
            return 0.;
    }

    // add the equilibrium emissivity of the population
    void addEquilibrium(Array& ev, double Teq) const
    {
        PlanckFunction B(Teq);
        int numLambda = _emlambdav.size();
        for (int ell = 0; ell != numLambda; ++ell)
        {
            ev[ell] += _emsigmaabsv[ell] * B(_emlambdav[ell]);
        }
    }

    // calculate the probabilities
    // Pv: the calculated probabilities (out)
    // ioff: the index offset in the temperature grid used for this calculation (out)
    // Am: scratch memory for the calculation (internal only)
    // Tmin/Tmax: temperature range in which to perform the calculation (in), and
    //            temperature range where the calculated probabilities are above a certain fraction of maximum (out)
    // Jv: the radiation field discretized on the input wavelength grid (in)
    void calcProbs(Array& Pv, int& ioff, Square<double>& Am, double& Tmin, double& Tmax, const Array& Jv) const
    {
        ioff = NR::locateClip(_grid->_Tv, Tmin);
        int NT = NR::locateClip(_grid->_Tv, Tmax) - ioff + 2;

        // copy/calculate the transition matrix coefficients
        Am.resize(NT);
        for (int f = 1; f < NT; f++)
        {
            const short* Kv = &_Km(f + ioff, ioff);
            const double* HRv = &_HRm(f + ioff, ioff);
            for (int i = 0; i < f; i++)
            {
                int k = Kv[i];
                Am(f, i) = k >= 0 ? HRv[i] * Jv[k] : 0.;
            }
        }
        for (int i = 1; i < NT; i++)
        {
            Am(i - 1, i) = _CRv[i + ioff];
        }

        // calculate the cumulative matrix coefficients, in place
        for (int f = NT - 2; f > 0; f--)
        {
            for (int i = 0; i < f; i++)
            {
                Am(f, i) += Am(f + 1, i);
            }
        }

        // calculate the probabilities
        Pv.resize(NT);
        Pv[0] = 1.;
        for (int i = 1; i < NT; i++)
        {
            double sum = 0.;
            for (int j = 0; j < i; j++) sum += Am(i, j) * Pv[j];
            Pv[i] = sum / Am(i - 1, i);

            // rescale if needed to keep infinities from happening
            if (Pv[i] > 1e10)
                for (int j = 0; j <= i; j++) Pv[j] /= Pv[i];
        }

        // normalize probabilities to unity
        Pv /= Pv.sum();

        // determine the temperature range where the probabability is above a given fraction of its maximum
        double frac = 1e-20 * Pv.max();
        int k;
        for (k = 0; k != NT - 2; k++)
            if (Pv[k] > frac) break;
        Tmin = _grid->_Tv[k + ioff];
        for (k = NT - 2; k != 1; k--)
            if (Pv[k] > frac) break;
        Tmax = _grid->_Tv[k + 1 + ioff];
    }

    // add the stochastic emissivity of the population (assumes that calcProbs has been called)
    // ev: the accumulated emissivity (in/out)
    // Tmin/Tmax: temperature range in which to add radiation (in)
    // Pv: the probabilities calculated previously by this calculator (in)
    // ioff: the index offset in the temperature grid used for that previous calculation (in)
    void addStochastic(Array& ev, double Tmin, double Tmax, const Array& Pv, int ioff) const
    {
        int imin = NR::locateClip(_grid->_Tv, Tmin);
        int imax = NR::locateClip(_grid->_Tv, Tmax);

        for (int i = imin; i <= imax; i++)
        {
            ev += _emsigmaabsv * _grid->_Bvv[i] * Pv[i - ioff];
        }
    }
};

////////////////////////////////////////////////////////////////////

// configuration constants
namespace
{
    // we construct three specialized temperature grids:
    // - a coarse grid for quickly determining the appropriate temperature range
    //   (with smaller bins in the lower temperature range, which sees action most often)
    // - a medium and a fine grid for performing the actual calculation (with equally-spaced bins)
    const double Tuppermax = 3000;     // the largest temperature taken into account
    const int NTA = 20;                // the number of grid points in the coarse grid
    const double ratioA = 500.;        // the ratio between the largest and smallest bins in the coarse grid
    const double widthB = 4.;          // the average width of bins in the medium grid, in K
    const double ratioB = 1.;          // the ratio between the largest and smallest bins in the medium grid
    const double widthC = 2.;          // the average width of bins in the fine grid, in K
    const double ratioC = 3.;          // the ratio between the largest and smallest bins in the fine grid
    const double deltaTmedium = 200.;  // the smallest temperature range for which the medium grid is used

    // considering the temperature range over which the probability is nonzero (or larger than a very small fraction);
    // we assume the population to be in equilibrium if one or both of the following conditions is true:
    // - the temperature range is smaller than a given delta-T (i.e. it resembles a delta function)
    // - the equilibrium temperature lies outside of the temperature range
    const double deltaTeq = 10.;  // the cutoff width of the temperature range
}

////////////////////////////////////////////////////////////////////

StochasticDustEmissionCalculator::~StochasticDustEmissionCalculator()
{
    for (auto calculator : _calculatorsA) delete calculator;
    for (auto calculator : _calculatorsB) delete calculator;
    for (auto calculator : _calculatorsC) delete calculator;
    delete _gridA;
    delete _gridB;
    delete _gridC;
}

////////////////////////////////////////////////////////////////////

void StochasticDustEmissionCalculator::precalculate(SimulationItem* item, const Array& lambdav, const Array& sigmaabsv,
                                                    string grainType, double bulkDensity, double meanMass,
                                                    const StoredTable<1>& enthalpy)
{
    auto config = item->find<Configuration>();

    // obtain the simulation's radiation field wavelength grid
    auto radiationFieldWLG = config->radiationFieldWLG();
    radiationFieldWLG->setup();

    // get the index of the bin being added
    int b = _calculatorsA.size();

    // perform initialization that needs to happen only once
    if (!b)
    {
        // copy the simulation's radiation field wavelength grid
        int n = radiationFieldWLG->numBins();
        _rflambdav.resize(n);
        _rfdlambdav.resize(n);
        for (int k = 0; k != n; ++k)
        {
            _rflambdav[k] = radiationFieldWLG->wavelength(k);
            _rfdlambdav[k] = radiationFieldWLG->effectiveWidth(k);
        }

        // if requested by the configuration, precalculate the CMB source term; otherwise the array remains empty
        if (config->includeHeatingByCMB())
        {
            PlanckFunction B(Constants::Tcmb() * (1. + config->redshift()));
            _Bcmbv.resize(n);
            for (int k = 0; k != n; ++k) _Bcmbv[k] = B(_rflambdav[k]);
        }

        // copy the simulation's dust emission wavelength grid
        auto dustEmissionWLG = config->dustEmissionWLG();
        dustEmissionWLG->setup();
        _emlambdav = dustEmissionWLG->extlambdav();

        // build the three temperature grids
        _gridA = new SDE_TemperatureGrid(_emlambdav, 2., Tuppermax, NTA, ratioA);
        _gridB = new SDE_TemperatureGrid(_emlambdav, 2., Tuppermax, static_cast<int>(Tuppermax / widthB), ratioB);
        _gridC = new SDE_TemperatureGrid(_emlambdav, 2., Tuppermax, static_cast<int>(Tuppermax / widthC), ratioC);
    }

    // build the three calculators for this bin and add them to the corresponding lists
    _calculatorsA.push_back(new SDE_Calculator(item, _gridA, radiationFieldWLG, _rflambdav, _rfdlambdav, _emlambdav,
                                               lambdav, sigmaabsv, bulkDensity, meanMass, enthalpy));
    _calculatorsB.push_back(new SDE_Calculator(item, _gridB, radiationFieldWLG, _rflambdav, _rfdlambdav, _emlambdav,
                                               lambdav, sigmaabsv, bulkDensity, meanMass, enthalpy));
    _calculatorsC.push_back(new SDE_Calculator(item, _gridC, radiationFieldWLG, _rflambdav, _rfdlambdav, _emlambdav,
                                               lambdav, sigmaabsv, bulkDensity, meanMass, enthalpy));

    // remember some other properties for this bin
    _meanMasses.push_back(meanMass);
    _grainTypes.push_back(grainType);
    _maxEnthalpyTemps.push_back(enthalpy.axisRange<0>().max());
}

////////////////////////////////////////////////////////////////////

size_t StochasticDustEmissionCalculator::allocatedBytes() const
{
    size_t allocatedBytes = 0;
    allocatedBytes += _rflambdav.size() * sizeof(_rflambdav[0]);
    allocatedBytes += _rfdlambdav.size() * sizeof(_rfdlambdav[0]);
    allocatedBytes += _Bcmbv.size() * sizeof(_Bcmbv[0]);
    allocatedBytes += _emlambdav.size() * sizeof(_emlambdav[0]);

    if (_gridA) allocatedBytes += _gridA->allocatedBytes();
    if (_gridB) allocatedBytes += _gridB->allocatedBytes();
    if (_gridC) allocatedBytes += _gridC->allocatedBytes();
    for (auto calculator : _calculatorsA) allocatedBytes += calculator->allocatedBytes();
    for (auto calculator : _calculatorsB) allocatedBytes += calculator->allocatedBytes();
    for (auto calculator : _calculatorsC) allocatedBytes += calculator->allocatedBytes();

    allocatedBytes += _meanMasses.size() * sizeof(_meanMasses[0]);
    allocatedBytes += _grainTypes.size() * sizeof(_grainTypes[0]);
    allocatedBytes += _maxEnthalpyTemps.size() * sizeof(_maxEnthalpyTemps[0]);
    return allocatedBytes;
}

////////////////////////////////////////////////////////////////////

Array StochasticDustEmissionCalculator::emissivity(const Array& Jv) const
{
    // if requested, create a local copy of the input radiation field that includes the CMB;
    // constructing a reference to either the input or this local copy avoids copying the input if there is no CMB
    Array Jcmbv;
    if (_Bcmbv.size()) Jcmbv = Jv + _Bcmbv;
    const Array& myJv = _Bcmbv.size() ? Jcmbv : Jv;

    // accumulate the emissivities in this array
    Array ev(_emlambdav.size());

    // this dictionary is updated as the loop over all bins in the mix proceeds;
    // for each type of grain composition, it keeps track of the grain mass above which
    // the representative grain is most certainly in equilibrium
    std::unordered_map<string, double> eqMass;

    // provide room for the probabilities calculated over each of the temperature grids
    Array Pv;
    Square<double> Am(static_cast<size_t>(ceil(Tuppermax)) + 1);

    // loop over all representative grains (size bins) in the dust mix
    int numBins = _calculatorsA.size();
    for (int b = 0; b != numBins; ++b)
    {
        // determine the equilibrium temperature for this bin using the calculator with a fine temperature grid
        double Teq = _calculatorsC[b]->equilibriumTemperature(myJv);

        // consider stochastic calculation only if the mean mass for this bin is below the cutoff mass
        string grainType = _grainTypes[b];
        double meanmass = _meanMasses[b];
        if (!eqMass.count(grainType) || meanmass < eqMass.at(grainType))
        {
            // calculate the probabilities over the coarse temperature grid
            double Tmin = 0;
            double Tmax = min(Tuppermax, _maxEnthalpyTemps[b]);

            int ioff = 0;
            _calculatorsA[b]->calcProbs(Pv, ioff, Am, Tmin, Tmax, myJv);

            // if the population might be stochastic...
            if (Tmax - Tmin > deltaTeq && Teq < Tmax)
            {
                // select the medium or fine temperature grid depending on the temperature range
                const SDE_Calculator* calculator = (Tmax - Tmin > deltaTmedium) ? _calculatorsB[b] : _calculatorsC[b];

                // calculate the probabilities over this grid, in the range determined by the coarse calculation
                calculator->calcProbs(Pv, ioff, Am, Tmin, Tmax, myJv);

                // if the population indeed is stochastic...
                if (Tmax - Tmin > deltaTeq && Teq < Tmax)
                {
                    // add the stochastic emissivity of this population to the running total
                    calculator->addStochastic(ev, Tmin, Tmax, Pv, ioff);
                    continue;
                }
            }

            // remember that all grains above this mass will be in equilibrium
            eqMass[grainType] = meanmass;
        }

        // otherwise, add the equilibrium emissivity of this population to the running total
        _calculatorsC[b]->addEquilibrium(ev, Teq);
    }
    return ev;
}

////////////////////////////////////////////////////////////////////
