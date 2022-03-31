/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MERIDIONALCUTFORM_HPP
#define MERIDIONALCUTFORM_HPP

#include "GenericForm.hpp"

//////////////////////////////////////////////////////////////////////

/** MeridionalCutForm represents a generic probe form. Refer to the ProbeFormBridge class for more
    information about probes and forms.

    This particular form outputs a text column file listing the quantity being probed at a number
    of positions along a meridian (i.e., one half of a great circle) through the spatial domain of
    the simulation. The half-circle is centered on the origin of the model coordinate system, lies
    in a meridional plane, and spans all inclinations (\f$0\leq\theta\leq\pi\f$). The radius
    \f$r\f$ and the azimuth \f$\varphi\f$ of the meridian can be configured, and the number of
    equidistant samples to be taken along the meridian can be specified as well.

    The output file contains a line for each sample along the meridian. The first column specifies
    the inclination \f$\theta\f$, and subsequent column(s) list the quantity being probed. */
class MeridionalCutForm : public GenericForm
{
    ITEM_CONCRETE(MeridionalCutForm, GenericForm, "a text column file with values along a meridian")
        ATTRIBUTE_TYPE_DISPLAYED_IF(MeridionalCutForm, "Level2")

        PROPERTY_INT(numSamples, "the number of samples along the meridian")
        ATTRIBUTE_MIN_VALUE(numSamples, "3")
        ATTRIBUTE_MAX_VALUE(numSamples, "100000")
        ATTRIBUTE_DEFAULT_VALUE(numSamples, "250")

        PROPERTY_DOUBLE(radius, "the radius of the circle containing the meridian")
        ATTRIBUTE_QUANTITY(radius, "length")
        ATTRIBUTE_MIN_VALUE(radius, "]0")

        PROPERTY_DOUBLE(azimuth, "the azimuth angle φ of the meridian")
        ATTRIBUTE_QUANTITY(azimuth, "posangle")
        ATTRIBUTE_MIN_VALUE(azimuth, "-360 deg")
        ATTRIBUTE_MAX_VALUE(azimuth, "360 deg")
        ATTRIBUTE_DEFAULT_VALUE(azimuth, "0 deg")
        ATTRIBUTE_DISPLAYED_IF(azimuth, "Dimension3")

    ITEM_END()

public:
    /** This function causes the form to output file(s) as described in the class header for the
        quantity being probed according to the information provided by the specified
        ProbeFormBridge instance. */
    void writeQuantity(const ProbeFormBridge* bridge) const override;
};

//////////////////////////////////////////////////////////////////////

#endif
