/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LISTGAUSSIANLINESSED_HPP
#define LISTGAUSSIANLINESSED_HPP

#include "GaussianLinesSED.hpp"

////////////////////////////////////////////////////////////////////

/** A ListGaussianLinesSED object represents a spectral energy distribution that consists of one
    or more Gaussian emission lines and that is fully specified inside the configuration file
    (i.e. without referring to an input file). It is intended for use in cases where there are just
    a few lines, but nothing keeps the user from specifying a long list. The lines may be listed in
    any order.

    The wavelengths are by default given in micron, but any wavelength, frequency or photon energy
    unit supported by SKIRT may be used instead. The velocity dispersions are given in velocity
    units. The luminosity values are given in luminosity units, but the scaling of the values is
    arbitrary because the %SED will be normalized after being loaded.

    See the description of the GaussianLinesSED class for more information on how the lines
    specified here are combined into a single tabulated %SED. */
class ListGaussianLinesSED : public GaussianLinesSED
{
    ITEM_CONCRETE(ListGaussianLinesSED, GaussianLinesSED,
                  "a multiple Gaussian-line SED specified inside the configuration file")

        PROPERTY_DOUBLE_LIST(wavelengths, "the line center wavelengths")
        ATTRIBUTE_QUANTITY(wavelengths, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelengths, "1 pm")
        ATTRIBUTE_MAX_VALUE(wavelengths, "1 m")

        PROPERTY_DOUBLE_LIST(dispersions, "the Gaussian velocity dispersion for each line")
        ATTRIBUTE_QUANTITY(dispersions, "velocity")
        ATTRIBUTE_MIN_VALUE(dispersions, "]0 km/s")
        ATTRIBUTE_MAX_VALUE(dispersions, "1000 km/s]")

        PROPERTY_DOUBLE_LIST(luminosities, "the relative luminosity of each line")
        ATTRIBUTE_QUANTITY(luminosities, "bolluminosity")
        ATTRIBUTE_MIN_VALUE(luminosities, "[0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function loads the configured wavelengths, dispersions and luminosities into the
        specified arrays. The GaussianLinesSED base class verifies that the three lists have a
        matching, nonzero length. */
    void getLineProperties(Array& lambdav, Array& dispersionv, Array& Lv) const override;
};

////////////////////////////////////////////////////////////////////

#endif
