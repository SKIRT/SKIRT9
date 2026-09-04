/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILEGAUSSIANLINESSED_HPP
#define FILEGAUSSIANLINESSED_HPP

#include "GaussianLinesSED.hpp"

////////////////////////////////////////////////////////////////////

/** A FileGaussianLinesSED object represents a spectral energy distribution that consists of one
    or more Gaussian emission lines and that is loaded from an input file. The floating point
    numbers in the first three columns of the text file specify, respectively, the center
    wavelength, the velocity dispersion, and the relative luminosity for each line. Any additional
    columns in the file are ignored. The lines may be listed in any order.

    The wavelengths are by default given in micron, the velocity dispersions in km/s, and the
    luminosities in solar units; all of these can be overridden by column header info in the file
    as usual for a TextInFile (this in particular means that the wavelength column may equally well
    list a frequency or a photon energy per line, using an appropriate unit). The scaling of the
    luminosity values is arbitrary because the %SED will be normalized after being loaded.

    See the description of the GaussianLinesSED class for more information on how the lines
    specified here are combined into a single tabulated %SED. */
class FileGaussianLinesSED : public GaussianLinesSED
{
    ITEM_CONCRETE(FileGaussianLinesSED, GaussianLinesSED, "a multiple Gaussian-line SED loaded from a text file")

        PROPERTY_STRING(filename, "the name of the file with the line definitions")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function loads the input file into the specified arrays. */
    void getLineProperties(Array& lambdav, Array& dispersionv, Array& Lv) const override;
};

////////////////////////////////////////////////////////////////////

#endif
