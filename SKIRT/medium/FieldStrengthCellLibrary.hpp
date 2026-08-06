/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FIELDSTRENGTHCELLLIBRARY_HPP
#define FIELDSTRENGTHCELLLIBRARY_HPP

#include "SpatialCellLibrary.hpp"

//////////////////////////////////////////////////////////////////////

/** The FieldStrengthCellLibrary class provides a library scheme for grouping spatial cells based
    on radiation field strength.

    The library contains a one-dimensional set of entries corresponding to different strengths of
    the stored radiation field, parameterized by the quantity \f[ U = \frac{ \int_0^\infty
    J_\lambda\, {\text{d}}\lambda }{ \int_0^\infty J_\lambda^{\text{MW}}\, {\text{d}}\lambda }, \f]
    where \f$J_\lambda^{\text{MW}}\f$ is the the local interstellar radiation field in the Milky
    Way according to Mathis et al. (1983). The mapping from spatial cells to library entries is
    built dynamically from binning the \f$N_{\text{cells}}\f$ values of \f$U\f$ as calculated from
    the medium system onto a one-dimensional grid with a user-configurable number of field strength
    bins. */
class FieldStrengthCellLibrary : public SpatialCellLibrary
{
    ITEM_CONCRETE(FieldStrengthCellLibrary, SpatialCellLibrary,
                  "a library scheme for grouping spatial cells based on radiation field strength")
        ATTRIBUTE_TYPE_INSERT(FieldStrengthCellLibrary, "NonIdentitySpatialCellLibrary")

        PROPERTY_INT(numFieldStrengths, "the number of field strength bins")
        ATTRIBUTE_MIN_VALUE(numFieldStrengths, "10")
        ATTRIBUTE_MAX_VALUE(numFieldStrengths, "10000000")
        ATTRIBUTE_DEFAULT_VALUE(numFieldStrengths, "1000")

    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function returns the number of entries in the library. In this class the function
        returns the user-configured number of field strength bins. */
    int numEntries() const override;

    /** This function returns a vector \em nv with length \f$N_{\text{cells}}\f$ that maps each
        cell index \f$m\f$ to the corresponding library entry index \f$n_m\f$. In this class the
        function loops over all spatial cells and calculates the strength \f$U_m\f$ of the
        radiation field for each of them. Based on these values, a logarithmic grid in \f$U\f$
        embracing these values is established, \f[ U_{(n)} = U_{\text{min}} \left( \frac{
        U_{\text{max}} }{ U_{\text{min}} } \right)^{n/N_U} \qquad n=0,\ldots,N_U, \f] where
        \f$U_{\text{min}}\f$ and \f$U_{\text{max}}\f$ represent the smallest and largest values of
        the field strength found among all cells. Then the function determines for each cell
        \f$m\f$ the corresponding library entry \f$n\f$. */
    vector<int> mapping(const Array& bv) const override;
};

////////////////////////////////////////////////////////////////////

#endif
