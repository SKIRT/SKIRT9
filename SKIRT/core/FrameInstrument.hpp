/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FRAMEINSTRUMENT_HPP
#define FRAMEINSTRUMENT_HPP

#include "DistantInstrument.hpp"

////////////////////////////////////////////////////////////////////

/** A FrameInstrument object represents a distant instrument that records the surface brightness in
    every pixel of a given frame for each wavelength, and outputs an IFU data cube in a FITS file.

    The instrument allows configuring the field of view and number of pixels in both directions of
    the observer plane. Photon packets arriving from a point that parallel projects outside of the
    field of view are ignored. */
class FrameInstrument : public DistantInstrument
{
    ITEM_CONCRETE(FrameInstrument, DistantInstrument,
                  "a distant instrument that outputs the surface brightness in every pixel as a data cube")

        PROPERTY_DOUBLE(fieldOfViewX, "the total field of view in the horizontal direction")
        ATTRIBUTE_QUANTITY(fieldOfViewX, "length")
        ATTRIBUTE_MIN_VALUE(fieldOfViewX, "]0")

        PROPERTY_INT(numPixelsX, "the number of pixels in the horizontal direction")
        ATTRIBUTE_MIN_VALUE(numPixelsX, "1")
        ATTRIBUTE_MAX_VALUE(numPixelsX, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numPixelsX, "250")

        PROPERTY_DOUBLE(centerX, "the center of the frame in the horizontal direction")
        ATTRIBUTE_QUANTITY(centerX, "length")
        ATTRIBUTE_DEFAULT_VALUE(centerX, "0")
        ATTRIBUTE_DISPLAYED_IF(centerX, "Level2")

        PROPERTY_DOUBLE(fieldOfViewY, "the total field of view in the vertical direction")
        ATTRIBUTE_QUANTITY(fieldOfViewY, "length")
        ATTRIBUTE_MIN_VALUE(fieldOfViewY, "]0")

        PROPERTY_INT(numPixelsY, "the number of pixels in the vertical direction")
        ATTRIBUTE_MIN_VALUE(numPixelsY, "1")
        ATTRIBUTE_MAX_VALUE(numPixelsY, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numPixelsY, "250")

        PROPERTY_DOUBLE(centerY, "the center of the frame in the vertical direction")
        ATTRIBUTE_QUANTITY(centerY, "length")
        ATTRIBUTE_DEFAULT_VALUE(centerY, "0")
        ATTRIBUTE_DISPLAYED_IF(centerY, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function configures the FluxRecorder instance associated with this instrument. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function simulates the detection of a photon packet by the instrument. It determines
        the projected position of the photon packet's last interaction site on the instrument frame
        and then calls the detect() function of the FluxRecorder instance associated with this
        instrument. */
    void detect(PhotonPacket* pp) override;

private:
    /** This private helper function returns the index of the spatial pixel on the detector that
        will be hit by a photon packet, or -1 if the photon packet does not hit the detector. Given
        the position \f${\boldsymbol{x}}=(x,y,z)\f$ of the last emission or scattering event of the
        photon packet, the direction \f${\boldsymbol{k}}_{\text{obs}} = (\theta,\varphi)\f$ towards
        the observer, and the roll angle \f$\omega\f$ of the instrument, the impact coordinates
        \f$(x_{\text{p}},y_{\text{p}})\f$ are given by the following Euler-like transformation,
        where \f$z_{\text{p}}\f$ is ignored: \f[ \begin{bmatrix}x_{\text{p}} \\ y_{\text{p}} \\
        z_{\text{p}} \end{bmatrix} = \begin{bmatrix}\cos\omega & -\sin\omega & 0\\ \sin\omega &
        \cos\omega & 0\\ 0 & 0 & 1 \end{bmatrix} \begin{bmatrix} 0 & 1 & 0 \\ -1 & 0 & 0 \\ 0 & 0 &
        1 \end{bmatrix} \begin{bmatrix} \cos\theta & 0 & -\sin\theta \\ 0 & 1 & 0 \\ \sin\theta & 0
        & \cos\theta \end{bmatrix} \begin{bmatrix}\cos\varphi & -\sin(-\varphi) & 0 \\
        \sin(-\varphi) & \cos\varphi & 0 \\ 0 & 0 & 1 \end{bmatrix} \begin{bmatrix}x\\y\\z
        \end{bmatrix} \f] In other words, the originating position is rotated about the Z-axis over
        the azimuth angle \f$\varphi\f$ (with a minus sign because the observer is looking towards
        the center rather than along the specified direction), then rotated about the new Y-axis
        over the inclination angle \f$\theta\f$, and finally rotated about the new Z-axis over the
        roll angle \f$\omega\f$ reduced by 90 degrees (this constant transformation over -90
        degrees is represented above as a separate matrix). The 90 degree correction on the
        roll angle is introduced so that it would be more natural to specify this angle; in
        most cases it can be left to its default value of 0. Given these impact coordinates, the
        pixel indices \f$i\f$ and \f$j\f$ are determined as \f[ \begin{split} i &=
        \frac{{\text{floor}}(x_{\text{p}}-x_{\text{min}})}{\Delta x} \\ j &=
        \frac{{\text{floor}}(y_{\text{p}}-y_{\text{max}})}{\Delta y} \end{split} \f] where
        \f${\text{floor}}(z)\f$ is an operator that returns the largest integer that is not greater
        than \f$y\f$. The spatial pixel number \f$l\f$ is then determined as \f$l=i+j\,N_x\f$,
        asuming \f$i\f$ and \f$j\f$ are indeed within the detector range. */
    int pixelOnDetector(const PhotonPacket* pp) const;

    //======================== Data Members ========================

private:
    // data members derived from the discoverable properties during setup, used in pixelOnDetector()
    double _costheta{0};
    double _sintheta{0};
    double _cosphi{0};
    double _sinphi{0};
    double _cosomega{0};
    double _sinomega{0};
    int _Nxp{0};
    int _Nyp{0};
    double _xpmin{0};
    double _xpsiz{0};
    double _ypmin{0};
    double _ypsiz{0};
};

////////////////////////////////////////////////////////////////////

#endif
