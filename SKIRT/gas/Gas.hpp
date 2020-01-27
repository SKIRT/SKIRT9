/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */
#ifndef GAS_HPP
#define GAS_HPP

#include "Array.hpp"
class WavelengthGrid;
namespace GasModule
{
    class GasInterface;
    class GasState;
}

#include <vector>

/** The stuff directly interfacing with the gas module should go into one or multiple classes of
    this subproject, so that the code can easily be built without the gas module. */
class Gas
{
public:
    Gas();
    ~Gas();
    void setup(const Array& frequencyv);
    void allocateGasStates(size_t num);
    void updateGasState(int m, const Array& meanIntensity);
    double gasTemperature(int m) const;
private:
    std::vector<GasModule::GasState> _statev;
    std::unique_ptr<GasModule::GasInterface> _gi;
};

#endif
