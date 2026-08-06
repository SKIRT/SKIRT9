/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TEMPERATUREWAVELENGTHCELLLIBRARY_HPP
#define TEMPERATUREWAVELENGTHCELLLIBRARY_HPP

#include "SpatialCellLibrary.hpp"

//////////////////////////////////////////////////////////////////////

/** The TemperatureWavelengthCellLibrary class provides a library scheme for grouping spatial cells
    based on the indicative dust temperature and the indicative dust wavelength of the stored
    radiation field in each cell. The library contains a two-dimensional table of entries
    corresponding to different values for each of these properties. The user can configure the
    number of bins in each of the two dimensions.

    Note: the current implementation of this class takes into account only the dust media; media
    with other material types are ignored.

    The indicative dust temperature for a cell is obtained by averaging the LTE equilibrium
    temperatures for the various dust mixes present in the cell. For a definition and more
    information about the indicative dust temperature, refer to the
    MediumSystem::indicativeDustTemperature() function.

    The indicative dust wavelength for a cell \f$m\f$ is defined as \f[ {\bar{\lambda}}_m =
    \frac{\int_0^\infty k_{m,\lambda}^{\text{abs}}\, J_{m,\lambda}\, \lambda\, {\text{d}}\lambda }{
    \int_0^\infty k_{m,\lambda}^{\text{abs}}\, J_{m,\lambda}\, {\text{d}}\lambda} \f] where
    \f$k_{m,\lambda}^{\text{abs}}\f$ is the dust absorption opacity in the cell, \f[
    k_{m,\lambda}^{\text{abs}} = \sum_h \kappa_{m,h,\lambda}^{\text{abs}} \, \rho_{m,h}\f] with the
    sum running over dusty media only. */
class TemperatureWavelengthCellLibrary : public SpatialCellLibrary
{
    ITEM_CONCRETE(TemperatureWavelengthCellLibrary, SpatialCellLibrary,
                  "a library scheme for grouping spatial cells based on indicative temperature and wavelength")
        ATTRIBUTE_TYPE_INSERT(TemperatureWavelengthCellLibrary, "NonIdentitySpatialCellLibrary")

        PROPERTY_INT(numTemperatures, "the number of temperature bins")
        ATTRIBUTE_MIN_VALUE(numTemperatures, "5")
        ATTRIBUTE_MAX_VALUE(numTemperatures, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numTemperatures, "40")

        PROPERTY_INT(numWavelengths, "the number of wavelength bins")
        ATTRIBUTE_MIN_VALUE(numWavelengths, "5")
        ATTRIBUTE_MAX_VALUE(numWavelengths, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numWavelengths, "25")

    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function returns the number of entries in the library. In this class the function
        returns the user-configured number of bins, i.e. the product of the number of temperature
        bins and the number of wavelength bins. */
    int numEntries() const override;

    /** This function returns a vector \em nv with length \f$N_{\text{cells}}\f$ that maps each
        cell index \f$m\f$ to the corresponding library entry index \f$n_m\f$. In this class the
        function loops over all spatial cells and calculates the indicative dust temperature and
        the indicative dust wavelength of the stored radiation field for each of them. Based on
        these values, a two-dimensional grid is established such that it fits all the measured
        values. The temperature grid points \f${\bar{T}}_{(i)}\f$ are distributed linearly, i.e.
        \f[ {\bar{T}}_{(i)} = {\bar{T}}_{\text{min}} + \frac{i}{N_{\bar{T}}}\,
        ({\bar{T}}_{\text{max}} - {\bar{T}}_{\text{min}}) \qquad i=0,\ldots,N_{\bar{T}} \f] where
        \f${\bar{T}}_{\text{min}}\f$ and \f${\bar{T}}_{\text{max}}\f$ represent the smallest and
        largest values of the indicative dust temperature found among all spatial cells. The
        wavelength grid points have logarithmic distribution, \f[ {\bar{\lambda}}_{(j)} =
        {\bar{\lambda}}_{\text{min}} \left( \frac{ {\bar{\lambda}}_{\text{max}} }{
        {\bar{\lambda}}_{\text{min}} } \right)^{j/N_{\bar{\lambda}}} \qquad
        j=0,\ldots,N_{\bar{\lambda}}. \f] The function then calculates for each cell \f$m\f$ its
        library entry \f$n \equiv (i,j)\f$. */
    vector<int> mapping(const Array& bv) const override;
};

////////////////////////////////////////////////////////////////////

#endif
