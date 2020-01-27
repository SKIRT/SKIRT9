#include "Gas.hpp"
#include "Array.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"

#ifdef BUILD_WITH_GAS
#    include "GasInterface.hpp"
#    include <iostream>
#endif

Gas::Gas()
{
#ifndef BUILD_WITH_GAS
    FATALERROR("SKIRT was built using BUILD_WITH_GAS=OFF, gas is not supported");
#endif
}

Gas::~Gas() = default;

void Gas::setup(const Array& frequencyv)
{
#ifdef BUILD_WITH_GAS
    _gi = std::make_unique<GasModule::GasInterface>(frequencyv, frequencyv, frequencyv);
#else
    (void)frequencyv;
#endif
}

void Gas::allocateGasStates(size_t num)
{
#ifdef BUILD_WITH_GAS
    _statev.resize(num);
#else
    (void)num;
#endif
}

void Gas::updateGasState(int m, const Array& meanIntensityv)
{
#ifdef BUILD_WITH_GAS
    const Array& iFrequencyv = _gi->iFrequencyv();
    Array jnu(meanIntensityv.size());

    if (iFrequencyv.size() != meanIntensityv.size())
        throw FATALERROR("Something went wrong with the wavelength/frequency grids");

    int countzeros = 0;
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

        if (jnu[i] <= 0)
        {
            jnu[i] = 1e-99;
            countzeros++;
        }

        // std::cout << "freq " << nu << " wav " << Constants::c() / nu << " jnu(erg cm-2) " << jnu[i]
        //           << " jlambda(J m-2) " << meanIntensityv[meanIntensityv.size() - 1 - i] << '\n';
    }
    std::cout << countzeros << " zeros in cell " << m << '\n';

    // No grains for now
    GasModule::GrainInterface gr;
    _gi->updateGasState(_statev[m], 1000., jnu, gr);
#else
    (void)m;
    (void)meanIntensityv;
#endif
}

double Gas::gasTemperature(int m) const
{
#ifdef BUILD_WITH_GAS
    return _statev[m].temperature();
#else
    (void)m;
    return 0;
#endif
}
