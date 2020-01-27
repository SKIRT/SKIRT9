/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */
#ifndef GAS_HPP
#define GAS_HPP

#include "Array.hpp"

/** The stuff directly interfacing with the gas module should go into one or multiple classes of
    this subproject, so that the code can easily be built without the gas module. */
class Gas
{
public:
    static void initialize(const Array& frequencyv);
    static void finalize();
    static void allocateGasStates(size_t num);
    static void updateGasState(int m, const Array& meanIntensity);
    static double gasTemperature(int m);
};

#endif
