/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StochasticDustEmissivity.hpp"
#include "ArrayTable.hpp"
#include "Constants.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "MultiGrainDustMix.hpp"
#include "NR.hpp"
#include "PlanckFunction.hpp"

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
        Square(size_t n) : _n(n), _v(new T[n*n]) { }
        ~Square() { delete[] _v; }

        // sets logical size, which must not be larger than maximum size set in constructor (is not checked)
        // does not clear values; does not resize underlying memory
        void resize(size_t n) { _n = n; }

        // access to values  (const version currently not needed)
        T& operator()(size_t i, size_t j) { return _v[i*_n+j]; }
    };

    // square matrix with only items below the diagonal (i>j)
    template<typename T> class Triangle
    {
    private:
        T* _v;
        static size_t offset(size_t i) { return ((i-1)*i)>>1; }

    public:
        // constructor sets size (can't be changed)
        Triangle(size_t n) : _v(new T[offset(n)]) { }
        ~Triangle() { delete[] _v; }

        // access to values; must have i>j (is not checked)
        const T& operator()(size_t i, size_t j) const { return _v[offset(i)+j]; }
        T& operator()(size_t i, size_t j) { return _v[offset(i)+j]; }
    };
}

////////////////////////////////////////////////////////////////////

// helper class to construct and store a particular temperature grid and the black-body emissivity spectrum
// discretized on this temperature grid and on the output wavelength grid;
// all members are public for ease of use in SDE_Calculator
class SDE_Grid
{
public:
    // output wavelength grid (indexed on ell)
    int _NoutLambda;
    Array _outLambdav;

    // temperature grid (indexed on i)
    int _NT;
    Array _Tv;

    // black-body radiation (indexed on i and ell)
    ArrayTable<2> _Bvv;

    // constructor
    SDE_Grid(const WavelengthGrid* outputWLG, double Tmin, double Tmax, int NT, double ratio)
    {
        // copy output wavelength grid
        _NoutLambda = outputWLG->numBins();
        _outLambdav.resize(_NoutLambda);
        for (int ell=0; ell!=_NoutLambda; ++ell) _outLambdav[ell] = outputWLG->wavelength(ell);

        // build temperature grid (linear if ratio==1)
        _NT = NT;
        NR::buildPowerLawGrid(_Tv, Tmin, Tmax, NT-1, ratio);

        // pre-calculate black-body radiation
        _Bvv.resize(NT,_NoutLambda);
        for (int i=0; i!=NT; ++i)
        {
            PlanckFunction B(_Tv[i]);
            for (int ell=0; ell!=_NoutLambda; ++ell) _Bvv(i,ell) = B(_outLambdav[ell]);
        }
    }
};

////////////////////////////////////////////////////////////////////

// helper class to calculate the temperature probability distribution on a particular fixed temperature grid
// for a representative grain embedded in a given radiation field (discretized on the input wavelength grid)
// with support for pre-calculating and storing relevant data relating to the grain/temperature grid combination
class SDE_Calculator
{
private:
    const SDE_Grid* _grid;      // temperature grid with corresponding black body discretization
    Triangle<double> _HRm;      // heating rates (indexed on f,i)
    Triangle<short> _Km;        // radiation field wavelength index k (indexed on f,i)
    Array _CRv;                 // cooling rates (indexed on i)
    Array _sigmaabsv;           // cross sections on output wavelength grid (indexed on ell)

public:
    SDE_Calculator(const SDE_Grid* grid, const WavelengthGrid* inputWLG, const MultiGrainDustMix* mgmix, int b)
        : _grid(grid), _HRm(grid->_NT), _Km(grid->_NT), _CRv(grid->_NT), _sigmaabsv(grid->_outLambdav)
    {
        // get shortcuts to the temperature grid
        const Array& Tv = grid->_Tv;
        int NT = grid->_NT;

        // calculate the enthalpy of a single dust grain in this population across the temperature grid
        Array Hv(NT);
        for (int i=0; i<NT; i++) Hv[i] = mgmix->binEnthalpy(b,Tv[i]);

        // calculate the enthalpy bin widths
        Array dHv(NT);
        dHv[0] = Hv[1]-Hv[0];
        for (int i=1; i<NT-1; i++)
        {
            double Tmin = (Tv[i-1]+Tv[i])/2.;
            double Tmax = (Tv[i+1]+Tv[i])/2.;
            dHv[i] = mgmix->binEnthalpy(b,Tmax) - mgmix->binEnthalpy(b,Tmin);
        }
        dHv[NT-1] = Hv[NT-1]-Hv[NT-2];

        // calculate the heating rates, barring the dependency on the radiation field
        double hc = Constants::h() * Constants::c();
        for (int f=1; f<NT; f++)
        {
            for (int i=0; i<f; i++)
            {
                double Hdiff = Hv[f] - Hv[i];
                double lambda = hc / Hdiff;
                int k = inputWLG->bin(lambda);
                if (k>=0) _HRm(f,i) = hc * mgmix->sectionAbs(lambda) * dHv[f] / (Hdiff*Hdiff*Hdiff);
                _Km(f,i) = k;
            }
        }

        // calculate the cooling rates
        const Array& lambdav = mgmix->finelambdav();
        const Array& sigmaabsv = mgmix->finesigmaabsv();
        int Nlambda = lambdav.size();
        for (int i=1; i<NT; i++)
        {
            PlanckFunction B(Tv[i]);
            double planckabs = 0.;
            for (int ell = 1; ell<Nlambda; ell++) // skip the first wavelength so we can determine a bin width
            {
                double lambda = lambdav[ell];
                double dlambda = lambdav[ell] - lambdav[ell-1];
                planckabs += sigmaabsv[ell] * B(lambda) * dlambda;
            }
            double Hdiff = Hv[i] - Hv[i-1];
            _CRv[i] = planckabs / Hdiff;
        }

        // obtain the cross sections on the output wavelength grid
        for (int ell=0; ell!=grid->_NoutLambda; ++ell)
        {
            _sigmaabsv[ell] = mgmix->binSectionAbs(b, grid->_outLambdav[ell]);
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
        for (int f=1; f<NT; f++)
        {
            const short* Kv = &_Km(f+ioff,ioff);
            const double* HRv = &_HRm(f+ioff,ioff);
            for (int i=0; i<f; i++)
            {
                int k = Kv[i];
                Am(f,i) = k>=0 ? HRv[i] * Jv[k] : 0.;
            }
        }
        for (int i=1; i<NT; i++)
        {
            Am(i-1,i) = _CRv[i+ioff];
        }

        // calculate the cumulative matrix coefficients, in place
        for (int f=NT-2; f>0; f--)
        {
            for (int i=0; i<f; i++)
            {
                Am(f,i) += Am(f+1,i);
            }
        }

        // calculate the probabilities
        Pv.resize(NT);
        Pv[0] = 1.;
        for (int i=1; i<NT; i++)
        {
            double sum = 0.;
            for (int j=0; j<i; j++) sum += Am(i,j) * Pv[j];
            Pv[i] = sum / Am(i-1,i);

            // rescale if needed to keep infinities from happening
            if (Pv[i] > 1e10) for (int j=0; j<=i; j++) Pv[j]/=Pv[i];
        }

        // normalize probabilities to unity
        Pv /= Pv.sum();

        // determine the temperature range where the probabability is above a given fraction of its maximum
        double frac = 1e-20 * Pv.max();
        int k;
        for (k=0; k!=NT-2; k++) if (Pv[k]>frac) break;
        Tmin = _grid->_Tv[k+ioff];
        for (k=NT-2; k!=1; k--) if (Pv[k]>frac) break;
        Tmax = _grid->_Tv[k+1+ioff];
    }

    // add the transient emissivity of the population
    // ev: the accumulated emissivity (in/out)
    // Tmin/Tmax: temperature range in which to add radiation (in)
    // Pv: the probabilities calculated previously by this calculator (in)
    // ioff: the index offset in the temperature grid used for that previous calculation (in)
    void addTransient(Array& ev, double Tmin, double Tmax, const Array& Pv, int ioff) const
    {
        int imin = NR::locateClip(_grid->_Tv, Tmin);
        int imax = NR::locateClip(_grid->_Tv, Tmax);

        for (int i=imin; i<=imax; i++)
        {
            ev += _sigmaabsv * _grid->_Bvv[i] * Pv[i-ioff];
        }
    }

    // add the equilibrium emissivity of the population
    void addEquilibrium(Array& ev, double Teq) const
    {
        PlanckFunction B(Teq);
        for (int ell=0; ell!=_grid->_NoutLambda; ++ell)
        {
            ev[ell] += _sigmaabsv[ell] * B(_grid->_outLambdav[ell]);
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
    const double Tuppermax = 3000;  // the largest temperature taken into account
    const int NTA = 20;             // the number of grid points in the coarse grid
    const double ratioA = 500.;     // the ratio between the largest and smallest bins in the coarse grid
    const double widthB = 4.;       // the average width of bins in the medium grid, in K
    const double ratioB = 1.;       // the ratio between the largest and smallest bins in the medium grid
    const double widthC = 2.;       // the average width of bins in the fine grid, in K
    const double ratioC = 3.;       // the ratio between the largest and smallest bins in the fine grid
    const double deltaTmedium = 200.;  // the smallest temperature range for which the medium grid is used

    // considering the temperature range over which the probability is nonzero (or larger than a very small fraction);
    // we assume the population to be in equilibrium if one or both of the following conditions is true:
    // - the temperature range is smaller than a given delta-T (i.e. it resembles a delta function)
    // - the equilibrium temperature lies outside of the temperature range
    const double deltaTeq = 10.;    // the cutoff width of the temperature range
}

////////////////////////////////////////////////////////////////////

StochasticDustEmissivity::~StochasticDustEmissivity()
{
    for (const auto& pairKeyValue : _calculatorsA) delete pairKeyValue.second;
    for (const auto& pairKeyValue : _calculatorsB) delete pairKeyValue.second;
    for (const auto& pairKeyValue : _calculatorsC) delete pairKeyValue.second;
    for (auto grid : _grids) delete grid;
}

////////////////////////////////////////////////////////////////////

void StochasticDustEmissivity::setupSelfBefore()
{
    DustEmissivity::setupSelfBefore();

    // get the input and output wavelength grids
    auto inputWLG = find<Configuration>()->radiationFieldWLG();
    auto outputWLG = find<Configuration>()->dustEmissionWLG();
    inputWLG->setup();
    outputWLG->setup();
    _numOutWavelengths = outputWLG->numBins();

    // construct the appropriate grids and calculators for each MultiGrainDustMix dust mix
    find<Log>()->info("Precalculating cached values for transient dust emissivity computations...");
    for (auto medium : find<MediumSystem>()->media())
    {
        // get the representative mix for this medium (TO DO: handle spatially varying dust mixes)
        auto mgmix = dynamic_cast<const MultiGrainDustMix*>(medium->mix());
        if (mgmix)
        {
            const_cast<MultiGrainDustMix*>(mgmix)->setup();  // TO DO: temporary hack

            // create grids
            double Tupper = min(Tuppermax, mgmix->maxEnthalpyTemperature());
            const SDE_Grid* gridA = new SDE_Grid(outputWLG, 2., Tupper, NTA, ratioA);
            const SDE_Grid* gridB = new SDE_Grid(outputWLG, 2., Tupper, static_cast<int>(Tupper/widthB), ratioB);
            const SDE_Grid* gridC = new SDE_Grid(outputWLG, 2., Tupper, static_cast<int>(Tupper/widthC), ratioC);
            _grids.push_back(gridA);
            _grids.push_back(gridB);
            _grids.push_back(gridC);

            // create calculators
            int numBins = mgmix->numBins();
            for (int b=0; b!=numBins; ++b)
            {
                _calculatorsA.emplace(std::make_pair(mgmix,b), new SDE_Calculator(gridA, inputWLG, mgmix, b));
                _calculatorsB.emplace(std::make_pair(mgmix,b), new SDE_Calculator(gridB, inputWLG, mgmix, b));
                _calculatorsC.emplace(std::make_pair(mgmix,b), new SDE_Calculator(gridC, inputWLG, mgmix, b));
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

Array StochasticDustEmissivity::emissivity(const MaterialMix* mix, const Array& Jv) const
{
    // accumulate the emissivities in this array
    Array ev(_numOutWavelengths);

    // we need the multi-grain dust mix interface; if this is not offered, return zero emissivity
    auto mgmix = dynamic_cast<const MultiGrainDustMix*>(mix);
    if (mgmix)
    {
        // this dictionary is updated as the loop over all bins in the mix proceeds;
        // for each type of grain composition, it keeps track of the grain mass above which
        // the representative grain is most certainly in equilibrium
        std::map<string,double> eqMass;

        // provide room for the probabilities calculated over each of the temperature grids
        Array Pv;
        Square<double> Am(static_cast<size_t>(ceil(Tuppermax))+1);

        // loop over all representative grains (size bins) in the dust mix
        int numBins = mgmix->numBins();
        for (int b=0; b!=numBins; ++b)
        {
            // determine the equilibrium temperature for this representative grain
            double Teq = mgmix->binEquilibriumTemperature(b, Jv);

            // get the coarse calculator for this representative grain
            const SDE_Calculator* calculatorA = _calculatorsA.at(std::make_pair(mgmix,b));

            // consider transient calculation only if the mean mass for this bin is below the cutoff mass
            string gcname = mgmix->binGrainType(b);
            double meanmass = mgmix->binMeanMass(b);
            if (!eqMass.count(gcname) || meanmass < eqMass.at(gcname))
            {
                // calculate the probabilities over the coarse temperature grid
                double Tmin = 0;
                double Tmax = Tuppermax;
                int ioff = 0;
                calculatorA->calcProbs(Pv, ioff, Am, Tmin, Tmax, Jv);

                // if the population might be transient...
                if (Tmax-Tmin > deltaTeq && Teq < Tmax)
                {
                    // select the medium or fine temperature grid depending on the temperature range
                    const SDE_Calculator* calculator = (Tmax-Tmin > deltaTmedium)
                                                       ? _calculatorsB.at(std::make_pair(mgmix,b))
                                                       : _calculatorsC.at(std::make_pair(mgmix,b));

                    // calculate the probabilities over this grid, in the range determined by the coarse calculation
                    calculator->calcProbs(Pv, ioff, Am, Tmin, Tmax, Jv);

                    // if the population indeed is transient...
                    if (Tmax-Tmin > deltaTeq && Teq < Tmax)
                    {
                        // add the transient emissivity of this population to the running total
                        calculator->addTransient(ev, Tmin, Tmax, Pv, ioff);
                        continue;
                    }
                }

                // remember that all grains above this mass will be in equilibrium
                eqMass[gcname] = meanmass;
            }

            // otherwise, add the equilibrium emissivity of this population to the running total
            calculatorA->addEquilibrium(ev, Teq);
        }

        // convert emissivity from "per hydrogen atom" to "per unit mass"
        ev /= mgmix->mass();
    }
    return ev;
}

////////////////////////////////////////////////////////////////////
