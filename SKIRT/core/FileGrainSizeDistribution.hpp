/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILEGRAINSIZEDISTRIBUTION_HPP
#define FILEGRAINSIZEDISTRIBUTION_HPP

#include "Array.hpp"
#include "GrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** FileGrainSizeDistribution represents a grain size distribution that is loaded from a text
    column input file. The floating point numbers in the first two columns of the text file specify
    respectively the grain size \f$a\f$ and the corresponding size distribution value
    \f$\text{dnda} \propto \frac{\text{d}n_\text{D}}{\text{d}a}\f$. Where needed, the tabulated
    values loaded from the file are interpolated logarithmically on both axes. Outside of the
    specified grain size range, the size distribution value is considered to be zero.

    The grain sizes \f$a\f$ are by default given in micron and must be listed in increasing order.
    The \f$\text{dnda}(a)\f$ values have untis of inverse length and are by default given in
    1/micron. These units can be overridden by column header info in the file. Any required
    normalization can and should be applied externally to this class in the configuration of the
    grain population (see GrainSizeDistribution and GrainPopulation). */
class FileGrainSizeDistribution : public GrainSizeDistribution
{
    ITEM_CONCRETE(FileGrainSizeDistribution, GrainSizeDistribution,
                  "a dust grain size distribution loaded from a text file")
        ATTRIBUTE_TYPE_DISPLAYED_IF(FileGrainSizeDistribution, "Level2")

        PROPERTY_STRING(filename, "the name of the file with the dust grain size distribution")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function loads the input file as described in the class header. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the minimum grain size \f$a_\text{min}\f$, i.e. the first grain size
        in the table loaded from file. */
    double amin() const override;

    /** This function returns the maximum grain size \f$a_\text{max}\f$, i.e. the last grain size
        in the table loaded from file. */
    double amax() const override;

    /** This function returns the value of the distribution \f$\text{dnda} \propto
        \frac{\text{d}n_\text{D}}{\text{d}a}\f$ for a given grain size \f$a\f$, log-log
        interpolated from the table loaded from file. If \f$a<a_\text{min}\f$ or
        \f$a>a_\text{max}\f$ the function returns zero. */
    double dnda(double a) const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Array _av;     // grain sizes
    Array _dndav;  // number of grain values
};

////////////////////////////////////////////////////////////////////

#endif
