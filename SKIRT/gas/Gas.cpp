#include "Gas.hpp"
#include "Array.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"
#include "GrainInterface.hpp"
#include "NR.hpp"
#include "ProcessManager.hpp"
#include "StringUtils.hpp"
#include "Table.hpp"

#ifdef BUILD_WITH_GAS
#    include "GasInterface.hpp"
#    include <chrono>
#    include <iostream>
#    include <vector>
namespace
{
    // set during initialize
    // ---------------------

    Array _lambdav;                         // wavelengths given at initialization
    Array _olambdav;                        // wavelengths for determining the index in the opacity table
    Array _elambdav;                        // wavelengths for calculating the emission
    GasModule::GasInterface* _gi;           // instance of GasInterface
    std::vector<Gas::DustInfo> _dustinfov;  // information about the dust populations in the simulation

    // set per cell during updateGasState
    // ----------------------------------
    std::vector<GasModule::GasState> _statev;  // result of the equilibrium calculation for each cell
    Table<2> _opacityvv;                       // opacity(m, ell) for each cell m and wavelength ell

    // utility functions
    // -----------------

    // translate a SKIRT grain type into one of the built-in grain types
    GasModule::GrainTypeLabel stringToGrainTypeLabel(const string& populationGrainType)
    {
        if (StringUtils::contains(populationGrainType, "Silicate"))
            return GasModule::GrainTypeLabel::SIL;
        else if (StringUtils::contains(populationGrainType, "Graphite")
                 || StringUtils::contains(populationGrainType, "PAH"))
            return GasModule::GrainTypeLabel::CAR;
        else
            return GasModule::GrainTypeLabel::OTHER;
    }

    // dlambda = - (c / nu^2) dnu
    // dnu = - (c / lambda^2) dlambda
    //
    // jnu dnu = jnu (- c / lambda^2) dlambda = jlambda dlambda
    // --> jnu = lambda^2 / c * jlambda
    // and conversion is involution (its own inverse)

    // convert from a quantity per x to a quantity per (c x^-1)
    Array x_to_cxm1(const Array& xv, const Array& quantity_xv)
    {
        int numx = xv.size();
        Array quantity_cxm1v(numx);
        for (int ix = 0; ix < numx; ix++)
        {
            double x = xv[ix];
            int icxm1 = numx - 1 - ix;
            quantity_cxm1v[icxm1] = x * x / Constants::c() * quantity_xv[ix];
        }
        return quantity_cxm1v;
    }

    Array lambdaToNu(const Array& lambdav, const Array& quantityPerLambda)
    {
        return x_to_cxm1(lambdav, quantityPerLambda);
    }

    Array nuToLambda(const Array& nuv, const Array& quantityPerNu) { return x_to_cxm1(nuv, quantityPerNu); }

    /** Convert array of x to array of 1 / x, with the elements ordered the other way round */
    Array invertAndFlip(const Array& xv)
    {
        int size = xv.size();
        Array xv_inv_flip(size);
        for (int i = 0; i < size; i++) xv_inv_flip[i] = 1. / xv[size - 1 - i];
        return xv_inv_flip;
    }

    // Get an array of grain number densities [cm-3] corresponding to dustInfo i with mix number
    // density mixNumberDens [m-3] (n in MediumSystem)
    Array mixNumberDensToGrainDensityv(int i, double mixNumberDens)
    {
        return _dustinfov[i].numberDensRatiov * mixNumberDens * 1.e-6;
    }

    // thread locals for efficiency
    // ----------------------------

    // properly initialized and modified by the first call to setThreadLocalGrainDensities
    thread_local GasModule::GrainInterface t_gr;
    thread_local bool t_gr_is_ready{false};
    void setThreadLocalGrainDensities(const Array& mixNumberDensv)
    {
        // initialize when a thread meets this function for the first time (i.e. no populations are
        // present yet)
        if (!t_gr_is_ready)
        {
            for (size_t i = 0; i < _dustinfov.size(); i++)
            {
                // Just use 30 as the initial guess for the dust temperature, since SKIRT doesn't really
                // support calculating the dust temperature for individual sizes.
                Array temperaturev(30., _dustinfov[i].sizev.size());
                // Set the grain number densities using the number density of the mix (fictional H
                // density), and change unit from m-3 to cm-3
                Array densityv = mixNumberDensToGrainDensityv(i, mixNumberDensv[i]);
                if (false)
                {
                    std::cout << "grain size:";
                    for (double d : _dustinfov[i].sizev) std::cout << ' ' << d;
                    std::cout << "\ngrain dens:";
                    for (double d : densityv) std::cout << ' ' << d;
                    std::cout << '\n';
                }
                t_gr.addPopulation(stringToGrainTypeLabel(_dustinfov[i].grainType), _dustinfov[i].sizev, densityv,
                                   temperaturev, _gi->iFrequencyv(), _dustinfov[i].qabsvv);
                t_gr_is_ready = true;
            }
        }
        else
        {
            // simply change the number densities of the populations added in the block above
            for (size_t i = 0; i < _dustinfov.size(); i++)
                t_gr.changePopulationDensityv(i, mixNumberDensToGrainDensityv(i, mixNumberDensv[i]));
        }
    }
}
#endif

////////////////////////////////////////////////////////////////////

void Gas::initialize(const Array& lambdav, const std::vector<DustInfo>& dustinfov, const Array& emissionWLG)
{
#ifdef BUILD_WITH_GAS
    if (_gi) FATALERROR("Gas module should be initialized exactly once");

    _lambdav = lambdav;
    _elambdav = emissionWLG;
    _dustinfov = dustinfov;

    // Change the units of the dust properties from SI to cgs
    for (Gas::DustInfo& d : _dustinfov)
    {
        // Change size unit m to cm
        d.sizev *= 100.;
        // Flip the qabs arrays, because frequencies. This happens in-place.
        for (size_t b = 0; b < d.qabsvv.size(); b++) std::reverse(std::begin(d.qabsvv[b]), std::end(d.qabsvv[b]));
    }

    // Calculate the input radiation field / output opacity frequency grid
    Array iFrequencyv = Constants::c() * invertAndFlip(_lambdav);

    // Calculate the output emissivity frequency grid
    Array eFrequencyv = Constants::c() * invertAndFlip(_elambdav);

    // derive a wavelength grid that will be used for converting a wavelength to an index in the
    // opacity table (copied from DustMix)
    int numLambda = _lambdav.size();
    _olambdav.resize(numLambda);
    _olambdav[0] = lambdav[0];
    for (int ell = 1; ell != numLambda; ++ell) _olambdav[ell] = sqrt(lambdav[ell] * lambdav[ell - 1]);

    // Turn off error handling (otherwise, gas module can call abort)
    GasModule::GasInterface::errorHandlersOff();
    // Initialize the gas module
    _gi = new GasModule::GasInterface(iFrequencyv, iFrequencyv, eFrequencyv);
#else
    (void)lambdav;
    (void)dustinfov;
    (void)emissionWLG;
#endif
}

////////////////////////////////////////////////////////////////////

void Gas::finalize()
{
#ifdef BUILD_WITH_GAS
    delete _gi;
    _gi = nullptr;
#endif
}

////////////////////////////////////////////////////////////////////

void Gas::allocateGasStates(size_t num)
{
#ifdef BUILD_WITH_GAS
    _statev.resize(num);
    _opacityvv.resize(num, _gi->oFrequencyv().size());
#else
    (void)num;
#endif
}

////////////////////////////////////////////////////////////////////

bool Gas::hasGrainTypeSupport(const string& populationGrainType)
{
    return stringToGrainTypeLabel(populationGrainType) != GasModule::GrainTypeLabel::OTHER;
}

////////////////////////////////////////////////////////////////////

void Gas::updateGasState(int m, double n, const Array& meanIntensityv, const Array& mixNumberDensv)
{
#ifdef BUILD_WITH_GAS
    auto start = std::chrono::high_resolution_clock::now();
    const Array& iFrequencyv = _gi->iFrequencyv();

    if (iFrequencyv.size() != meanIntensityv.size())
        throw FATALERROR("Something went wrong with the wavelength/frequency grids");

    Array jnu = lambdaToNu(_lambdav, meanIntensityv);
    // unit conversion:
    // for gas module: erg s-1 cm-2 sr-1 Hz-1
    // for skirt     : J   s-1 m-2  sr-1 Hz-1
    //                 7   0   -4
    jnu *= 1.e3;

    size_t countzeros = 0;
    for (size_t i = 0; i < iFrequencyv.size(); i++)
        if (jnu[i] <= 0) countzeros++;

    bool verbose = !(m % 300);
    if (verbose && countzeros) std::cout << countzeros << " zeros in cell " << m << '\n';

    // prepare grain info for this cell
    setThreadLocalGrainDensities(mixNumberDensv);

    // calculate the equilibrium
    _gi->updateGasState(_statev[m], n * 1.e-6, jnu, t_gr);

    // calculate and store the opacity; the opacity table is indexed on wavelength, so we need to
    // flip the result around
    const Array& opacity_nu = _gi->opacity(_statev[m], true);
    for (size_t ell = 0; ell < opacity_nu.size(); ell++) _opacityvv(m, ell) = opacity_nu[opacity_nu.size() - 1 - ell];

    if (verbose)
    {
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "gas sample " << m << " n " << n * 1.e-6 << " time " << duration.count() << " ms.\n";
        std::cout << _gi->quickInfo(_statev[m], jnu) << '\n';
    }
#else
    (void)n;
    (void)m;
    (void)meanIntensityv;
    (void)mixNumberDensv;
#endif
}

////////////////////////////////////////////////////////////////////

void Gas::communicateResults()
{
#ifdef BUILD_WITH_GAS
    ProcessManager::sumToAll(_opacityvv.data());
#endif
}

////////////////////////////////////////////////////////////////////

void Gas::clearResults()
{
#ifdef BUILD_WITH_GAS
    _opacityvv.setToZero();
#endif
}

////////////////////////////////////////////////////////////////////

double Gas::gasTemperature(int m)
{
#ifdef BUILD_WITH_GAS
    return _statev[m].temperature();
#else
    (void)m;
    return 0;
#endif
}

////////////////////////////////////////////////////////////////////

double Gas::opacityAbs(double lambda, int m)
{
#ifdef BUILD_WITH_GAS
    return _opacityvv(m, indexForLambda(lambda));
#else
    (void)lambda;
    (void)m;
    return 0;
#endif
}

////////////////////////////////////////////////////////////////////

double Gas::opacityAbs(int ell, int m)
{
#ifdef BUILD_WITH_GAS
    return _opacityvv(m, ell);
#else
    (void)ell;
    (void)m;
    return 0;
#endif
}

////////////////////////////////////////////////////////////////////

int Gas::indexForLambda(double lambda)
{
#ifdef BUILD_WITH_GAS
    return NR::locateClip(_olambdav, lambda);
#else
    (void)lambda;
    return 0;
#endif
}

////////////////////////////////////////////////////////////////////

Array Gas::emissivity(int m)
{
#ifdef BUILD_WITH_GAS
    return nuToLambda(_gi->eFrequencyv(), _gi->emissivity(_statev[m], true));
#else
    (void)m;
    return Array;
#endif
}

////////////////////////////////////////////////////////////////////


