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

#ifdef BUILD_WITH_GAS
#    include <vector>
#endif

/** The stuff directly interfacing with the gas module should go into one or multiple classes of
    this subproject, so that the code can easily be built without the gas module. Things to think
    about:

    - Whether to store the gas states for every cell, or to pick and choose information from them
      (temperature / opacity) and store that info in an easier to manage data structure.

    - Updating gas states (calling the interface)

    - Exctracting information from gas states

    - Passing the wavelength grid at setup

    - Passing the grain properties when updating gas states

    - Multiprocessing

    - Should anything be integrated with the MaterialMix, or maybe the future 'constituent' concept?
*/
class Gas
{
public:
    Gas();
    ~Gas();
    void setup(const Array& frequencyv);
    void allocateGasStates(size_t num);
    void updateGasState(int m, const Array& meanIntensity);
    double gasTemperature(int m) const;
#ifdef BUILD_WITH_GAS
private:
    std::vector<GasModule::GasState> _statev;
    std::unique_ptr<GasModule::GasInterface> _gi;
#endif
};

#endif
