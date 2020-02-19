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
    struct DustInfo
    {
        string grainType;
        // Representative sizes. Will be used naively in the gas code (processes are calculated
        // separately for each size given here, and then summed). Integrating everything over the
        // grain size distribution is not doable for most processes in the gas code.
        Array sizev;
        // Number of grains 'per hydrogen atom' for each size
        Array numberDensRatiov;
        // Q_abs(a, nu), indexed on (size, frequency)
        std::vector<Array> qabsvv;
    };

    /** The wavelength grid passed here should be the grid used for meanIntensityv. */
    static void initialize(const Array& lambdav, const std::vector<DustInfo>&);

    /** Function which should be called exactly one, before the code finishes. After calling
        finalize(), all info about the gas will be lost. */
    static void finalize();
    static void allocateGasStates(size_t num);

    /** This function returns \c true if the grain type described by the given string (get it from
        MultiGrainDustMic::populationGrainType(c)), is supported by the gas module. */
    static bool hasGrainTypeSupport(const string& populationGrainType);

    /** This function prepares the arguments for \c GasInterface::updateGasState() and calls it.
        The \c m argument indicates which of the allocated gas states should be updated. The number
        density of the gas medium and the radiation field (mean intensity) should be given as
        arguments, as well as the number densities for the dust media. The dust number densities
        for each size will be calculated by multiplying the numbers in \c mixNumberDensv with the
        \c numberDensRatiov of each \c DustInfo. This argument should contain the hydrogen number
        density factor \c n for medium \c h of cell \c m. */
    static void updateGasState(int m, double n, const Array& meanIntensityv, const Array& mixNumberDensv);

    static double gasTemperature(int m);
};

#endif
