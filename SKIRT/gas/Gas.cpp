#include "Gas.hpp"
#include "Array.hpp"

#ifdef BUILD_WITH_GAS
#    include "GasInterface.hpp"
#    include <iostream>
#endif

Gas::Gas()
{
#ifdef BUILD_WITH_GAS
    Array frequencyv;
    GasModule::GasInterface gi(frequencyv, frequencyv, frequencyv);
    // Do something, will probably crash
    GasModule::GasState gs;
    GasModule::GrainInterface gr;
    std::cout << "Running";
    gi.updateGasState(gs, 1000., frequencyv, gr, nullptr);
#endif
}
