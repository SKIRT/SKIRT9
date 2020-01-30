#include "Gas.hpp"
#include "Array.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"

#ifdef BUILD_WITH_GAS
#    include "GasInterface.hpp"
#    include <iostream>
#    include <vector>
namespace
{
    GasModule::GasInterface* _gi;
    std::vector<GasModule::GasState> _statev;
    std::vector<Gas::DustInfo> _dustinfov;
}
#endif

void Gas::initialize(const Array& frequencyv, const std::vector<DustInfo>& dustinfov)
{
#ifdef BUILD_WITH_GAS
    if (_gi) FATALERROR("Gas module should be initialized exactly once");
    _gi = new GasModule::GasInterface(frequencyv, frequencyv, frequencyv);
    _dustinfov = dustinfov;
#else
    (void)frequencyv;
#endif
}

void Gas::finalize()
{
#ifdef BUILD_WITH_GAS
    delete _gi;
    _gi = nullptr;
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

void Gas::updateGasState(int m, const Array& meanIntensityv, const Array& mixNumberDensv)
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

    // Make grain interface object
    GasModule::GrainInterface gr;
    for (size_t i = 0; i < _dustinfov.size(); i++)
    {
        GasModule::GrainTypeLabel type;
        if (_dustinfov[i].type == 1)
        {
            type = GasModule::GrainTypeLabel::SIL;
        }
        else if (_dustinfov[i].type == 2)
        {
            type = GasModule::GrainTypeLabel::CAR;
        }
        else
        {
            type = GasModule::GrainTypeLabel::OTHER;
            FATALERROR("Unsupported grain type for gas");
        }

        Array temperaturev(_dustinfov[i].sizev.size());
        Array densityv = _dustinfov[i].numberDensRatiov * mixNumberDensv[i];

        gr.addPopulation(type, _dustinfov[i].sizev, densityv, temperaturev, _gi->iFrequencyv(), _dustinfov[i].qabsvv);
    }

    _gi->updateGasState(_statev[m], 1000., jnu, gr);
#else
    (void)m;
    (void)meanIntensityv;
#endif
}

double Gas::gasTemperature(int m)
{
#ifdef BUILD_WITH_GAS
    return _statev[m].temperature();
#else
    (void)m;
    return 0;
#endif
}
