/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CLEARDENSITYRECIPE_HPP
#define CLEARDENSITYRECIPE_HPP

#include "DynamicStateRecipe.hpp"
#include <atomic>
class Log;

////////////////////////////////////////////////////////////////////

/** ClearDensityRecipe is a dynamic medium state recipe mainly intended for testing purposes. The
    recipe sets the density to zero for all medium components in all spatial cells where the
    radiation field strength \f$U\f$ exceeds a given threshold value \f$U_\mathrm{cut}\f$
    configured by the user. The recipe converges when no additional cells have been cleared after
    an iteration.

    The radiation field strength \f$U\f$ in this recipe is defined as \f[ U = \frac{ \int_0^\infty
    J_\lambda\, {\text{d}}\lambda }{ \int_0^\infty J_\lambda^{\text{MW}}\, {\text{d}}\lambda }, \f]
    where \f$J_\lambda^{\text{MW}}\f$ is the local interstellar radiation field in the Milky Way
    according to Mathis et al. (1983). */
class ClearDensityRecipe : public DynamicStateRecipe
{
    ITEM_CONCRETE(ClearDensityRecipe, DynamicStateRecipe,
                  "a recipe that clears material in cells above a given radiation field strength")

        PROPERTY_DOUBLE(fieldStrengthThreshold, "the field strength above which material is cleared from a cell")
        ATTRIBUTE_MIN_VALUE(fieldStrengthThreshold, "[0")
        ATTRIBUTE_DEFAULT_VALUE(fieldStrengthThreshold, "1")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function is called before an update cycle begins. It initializes the number of cleared
        cells to zero and caches the radiation field wavelength grid. */
    void beginUpdate() override;

    /** This function is called repeatedly as part of the update cycle. If the density of the
        medium component in the cell under study is nonzero, it calculates the strength of the
        radiation field. If this value exceeds the configured threshold, the density is set to zero
        and the counter of cleared cells is atomically incremented. */
    void update(MaterialState* state, const Array& Jv) override;

    /** This function is called at the end of the update cyle. It issues a log message reporting
        the number and percentage of cells cleared during this update cyle and returns true if no
        cells have been cleared during this cycle or false otherwise. */
    bool endUpdate() override;

    //======================== Data Members =======================

private:
    Array _dlambdav;                     // wavelength bin widths for the radiation field wavelength grid
    std::atomic<size_t> _numCleared{0};  // number of cell/component combinations cleared during this update cycle
};

////////////////////////////////////////////////////////////////////

#endif
