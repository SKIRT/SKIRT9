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
        // Representative sizes. Will be used naively in the gas code (processes are calculated
        // separately for each size given here, and then summed). Integrating everything over the
        // grain size distribution is not doable for most processes in the gas code.
        Array sizev;
        // Multiplying these numbers with the total number density of the medium (h) this dust
        // component (c) originates from, should yield the number density for each representative
        // size (a).
        Array numberDensRatiov;
        // Q_abs(a, nu), indexed on (size, frequency)
        std::vector<Array> qabsvv;
    };

    /** The wavelength grid passed here should be the grid used for meanIntensityv. */
    static void initialize(const Array& frequencyv, const std::vector<DustInfo>&);

    /** Function which should be called exactly one, before the code finishes. After calling
        finalize(), all info about the gas will be lost. */
    static void finalize();
    static void allocateGasStates(size_t num);

    /** The dust number densities for each size (the info that was passed during initialize) will
        be rescaled using the total number density of the mix. The latter depends on the cell. */
    static void updateGasState(int m, const Array& meanIntensityv, const Array& mixNumberDensv);

    static double gasTemperature(int m);
};

#endif
