/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */
#ifndef GAS_HPP
#define GAS_HPP

#include "Array.hpp"
#include "Table.hpp"

/** The stuff directly interfacing with the gas module should go into one or multiple classes of
    this subproject, so that the code can easily be built without the gas module. */
class Gas
{
public:
    struct DustInfo
    {
        // change to proper enum later. Now use 1 for silicate, 2 for graphite
        int type;
        Array sizev;
        Array nPerMassUnitv;  // number density per mass density unit (passed to updateGasState), for each size bin
        std::vector<Array> qabsvv;
    };
    static void initialize(const Array& frequencyv, const std::vector<DustInfo>&);
    static void finalize();
    static void allocateGasStates(size_t num);
    // Pass dust mass density per component here
    static void updateGasState(int m, const Array& meanIntensity, const Array& dustMassDensv);
    static double gasTemperature(int m);
};

#endif
