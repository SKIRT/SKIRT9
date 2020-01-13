/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MarastonSED.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"
#include "MarastonSEDFamily.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

const SEDFamily* MarastonSED::getFamilyAndParameters(Array& parameters)
{
    // verify that the configured parameter values are in the valid portion of grid
    if ((_metallicity < 0.000894 || _metallicity > 0.0447) && (_age < 1e9 * Constants::year()))
        throw FATALERROR("For metallicities Z<0.000894 or Z>0.0447, the age should be at least 1 Gyr");

    // set the parameters using arbitrary scaling
    NR::assign(parameters, 1., _metallicity, _age);

    // construct the library of SED models
    return new MarastonSEDFamily(this, _imf == IMF::Kroupa ? MarastonSEDFamily::IMF::Kroupa
                                                           : MarastonSEDFamily::IMF::Salpeter);
}

////////////////////////////////////////////////////////////////////
