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
    // Single instance of GasInterface managed by initialize and finalize
    GasModule::GasInterface* _gi;

    // Stores some minimal info about the gas equilibrium for each cell. Can be used to
    // (re-)calculate the opacity and emissivity of the gas.
    std::vector<GasModule::GasState> _statev;

    // Information about the dust populations in the simulation. Set by initialize, should not be
    // modified afterwards.
    std::vector<Gas::DustInfo> _dustinfov;

    // Translate a SKIRT grain type into one of the built-in grain types
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

    // Wavelengths for the opacities
    Array _olambdav;

    // Store the opacity in each cell in one big table, for easy communication between processes.
    // Indexed on m, ell.
    Table<2> _opacityvv;
}
#endif

////////////////////////////////////////////////////////////////////

void Gas::initialize(const Array& lambdav, const std::vector<DustInfo>& dustinfov)
{
#ifdef BUILD_WITH_GAS
    if (_gi) FATALERROR("Gas module should be initialized exactly once");

    // Turn off error handling (otherwise, gas module can call abort)
    GasModule::GasInterface::errorHandlersOff();

    // Calculate the frequency grid
    int numFreq = lambdav.size();
    Array frequencyv(numFreq);
    for (int i = 0; i < numFreq; i++) frequencyv[i] = Constants::c() / lambdav[numFreq - 1 - i];

    // Initialize the gas module using this frequency grid
    _gi = new GasModule::GasInterface(frequencyv, frequencyv, frequencyv);

    // Copy the dust properties
    _dustinfov = dustinfov;

    // Change the units of the dust properties from SI to cgs
    for (Gas::DustInfo& d : _dustinfov)
    {
        // Change size unit m to cm
        d.sizev *= 100.;
        // Flip the qabs arrays, because frequencies. This happens in-place.
        for (size_t b = 0; b < d.qabsvv.size(); b++) std::reverse(std::begin(d.qabsvv[b]), std::end(d.qabsvv[b]));
    }

    // derive a wavelength grid that will be used for converting a wavelength to an index in the
    // opacity table. copied from DustMix
    _olambdav.resize(numFreq);
    _olambdav[0] = lambdav[0];
    for (int ell = 1; ell != numFreq; ++ell)
    {
        _olambdav[ell] = sqrt(lambdav[ell] * lambdav[ell - 1]);
    }
#else
    (void)lambdav;
    (void)dustinfov;
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
    Array jnu(meanIntensityv.size());

    if (iFrequencyv.size() != meanIntensityv.size())
        throw FATALERROR("Something went wrong with the wavelength/frequency grids");

    size_t countzeros = 0;
    for (size_t i = 0; i < iFrequencyv.size(); i++)
    {
        double nu = iFrequencyv[i];
        // j_nu = lambda * lambda * j_lambda / c = c2 nu-2 j_lambda c-1 = c nu-2 j_lambda
        jnu[i] = Constants::c() / nu / nu * meanIntensityv[meanIntensityv.size() - 1 - i];

        // unit conversion:
        // for gas module: erg s-1 cm-2 sr-1 Hz-1
        // for skirt     : J   s-1 m-2  sr-1 Hz-1
        //                 7   0   -4
        jnu[i] *= 1e3;
        if (jnu[i] <= 0) countzeros++;
    }
    bool verbose = !(m % 300);
    if (verbose && countzeros) std::cout << countzeros << " zeros in cell " << m << '\n';

    // Make grain interface object
    GasModule::GrainInterface gr;
    for (size_t i = 0; i < _dustinfov.size(); i++)
    {
        // Just use 30 as the initial guess for the dust temperature, since SKIRT doesn't really
        // support calculating the dust temperature for individual sizes.
        Array temperaturev(30., _dustinfov[i].sizev.size());
        // Set the grain number densities using the number density of the mix (fictional H
        // density), and change unit from m-3 to cm-3
        Array densityv = _dustinfov[i].numberDensRatiov * mixNumberDensv[i] * 1.e-6;
        if (false && verbose)
        {
            std::cout << "grain size:";
            for (double d : _dustinfov[i].sizev) std::cout << ' ' << d;
            std::cout << "\ngrain dens:";
            for (double d : densityv) std::cout << ' ' << d;
            std::cout << '\n';
        }
        gr.addPopulation(stringToGrainTypeLabel(_dustinfov[i].grainType), _dustinfov[i].sizev, densityv, temperaturev,
                         _gi->iFrequencyv(), _dustinfov[i].qabsvv);
    }
    _gi->updateGasState(_statev[m], n * 1.e-6, jnu, gr);

    // if (m == 0)
    // {
    //     // write out nu Jnu for first cell
    //     double c_um = 2.99792458e14;
    //     Array lambda_um = c_um / _gi->iFrequencyv();
    //     for (size_t i = 0; i < lambda_um.size(); i++)
    //         std::cout << lambda_um[i] << " " << _gi->iFrequencyv()[i] * jnu[i] << '\n';
    // }

    // calculate and store the opacity
    const Array& opacity_nu = _gi->opacity(_statev[m], true);
    // the opacity table is indexed on wavelength, so we need to flip the result around
    for (size_t ell = 0; ell < opacity_nu.size(); ell++) _opacityvv(m, ell) = opacity_nu[opacity_nu.size() - 1 - ell];

    if (verbose)
    {
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "gas sample " << m << " n " << n * 1.e-6 << " time "
                  << duration.count() << " ms.\n";
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
