/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CastelliKuruczSED.hpp"
#include "CastelliKuruczSEDFamily.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

const SEDFamily* CastelliKuruczSED::getFamilyAndParameters(Array& parameters)
{
    // cutoff values for temperature and gravity (see table in class documentation)
    Array Tv = {49000, 39000, 31000, 26000, 19000, 11750, 9000, 8250, 7500, 6000};
    Array gv = {5, 4.5, 4, 3.5, 3, 2.5, 2, 1.5, 1, 0.5};
    gv = pow(10., gv - 2.);

    // verify that the configured parameter values are in the valid portion of the grid
    size_t n = Tv.size();
    for (size_t i = 0; i != n; ++i)
    {
        if (_temperature > Tv[i] && _gravity < gv[i])
            throw FATALERROR("If the temperature is above " + StringUtils::toString(Tv[i])
                             + " K then the gravity must be above " + StringUtils::toString(gv[i]) + " km/s2");
    }

    // set the parameters using arbitrary scaling
    NR::assign(parameters, 1., _metallicity, _temperature, _gravity);

    // construct the library of SED models
    return new CastelliKuruczSEDFamily(this);
}

////////////////////////////////////////////////////////////////////
