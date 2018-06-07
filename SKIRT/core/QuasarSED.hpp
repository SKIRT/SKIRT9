/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef QUASARSED_HPP
#define QUASARSED_HPP

#include "ResourceSED.hpp"

////////////////////////////////////////////////////////////////////

/** QuasarSED represents a simple model for the spectral energy distribution of a typical quasar;
    see Stalevski et al. (2012, MNRAS, 420, 2756–2772) and Schartmann et al. (2005, A&A, 437,
    861–881). It is defined in the wavelength range between 0.001 \f$\mu\f$m and 1000 \f$\mu\f$m
    and is characterized by \f[ S_\lambda \propto \begin{cases} \; \lambda^{1/5} & \text{if
    $0.001~\mu{\text{m}}<\lambda<0.01~\mu{\text{m}}$} \\ \; \lambda^{-1} & \text{if
    $0.01~\mu{\text{m}}<\lambda<0.1~\mu{\text{m}}$} \\ \; \lambda^{-3/2} & \text{if
    $0.1~\mu{\text{m}}<\lambda<5~\mu{\text{m}}$} \\ \; \lambda^{-4} & \text{if
    $5~\mu{\text{m}}<\lambda<1000~\mu{\text{m}}$.} \end{cases} \f] */
class QuasarSED : public ResourceSED
{
    ITEM_CONCRETE(QuasarSED, ResourceSED, "the spectral energy distribution of a typical quasar")
    ITEM_END()

    //======================== Other Functions =======================

    /** This function returns the name of the stored table resource tabulating the %SED. */
    string resourceName() const override;
};

////////////////////////////////////////////////////////////////////

#endif
