/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILEPOLARIZEDPOINTSOURCE_HPP
#define FILEPOLARIZEDPOINTSOURCE_HPP

#include "LuminosityNormalization.hpp"
#include "Source.hpp"
#include "VelocityInterface.hpp"
class ContSED;

//////////////////////////////////////////////////////////////////////

/** FilePolarizedPointSource represents a primary source limited to a single point in space and
    emitting polarized radiation with an axi-symmetric angular dependence. The Stokes vector
    components of the emitted radiation, as a function of wavelength and inclination angle, are
    loaded from a user-provided file. The input file specifies these quantities with arbitrary
    normalization; the bolometric output is characterized by a LuminosityNormalization object
    configured in the ski file.

    The class offers properties to specify the position of the point source, the orientation of the
    symmetry axis of the angular dependence, and a "bulk" velocity. */
class FilePolarizedPointSource : public Source, public VelocityInterface
{
    ITEM_CONCRETE(FilePolarizedPointSource, Source, "a primary point source with a polarized spectrum read from file")
        ATTRIBUTE_TYPE_DISPLAYED_IF(FilePolarizedPointSource, "Level2")
        ATTRIBUTE_TYPE_INSERT(FilePolarizedPointSource, "Dimension2")

        PROPERTY_STRING(filename, "the name of the stored table file listing the Stokes vector components")

        PROPERTY_ITEM(normalization, LuminosityNormalization, "the type of luminosity normalization for the source")
        ATTRIBUTE_DEFAULT_VALUE(normalization, "IntegratedLuminosityNormalization")

        PROPERTY_DOUBLE(positionX, "the position of the point source, x component")
        ATTRIBUTE_QUANTITY(positionX, "length")
        ATTRIBUTE_DEFAULT_VALUE(positionX, "0")
        ATTRIBUTE_INSERT(positionX, "positionX:Dimension3")

        PROPERTY_DOUBLE(positionY, "the position of the point source, y component")
        ATTRIBUTE_QUANTITY(positionY, "length")
        ATTRIBUTE_DEFAULT_VALUE(positionY, "0")
        ATTRIBUTE_INSERT(positionY, "positionY:Dimension3")

        PROPERTY_DOUBLE(positionZ, "the position of the point source, z component")
        ATTRIBUTE_QUANTITY(positionZ, "length")
        ATTRIBUTE_DEFAULT_VALUE(positionZ, "0")

        PROPERTY_DOUBLE(symmetryX, "the direction of the positive symmetry axis, x component")
        ATTRIBUTE_QUANTITY(symmetryX, "length")
        ATTRIBUTE_DEFAULT_VALUE(symmetryX, "0")
        ATTRIBUTE_INSERT(symmetryX, "symmetryX:Dimension3")

        PROPERTY_DOUBLE(symmetryY, "the direction of the positive symmetry axis, y component")
        ATTRIBUTE_QUANTITY(symmetryY, "length")
        ATTRIBUTE_DEFAULT_VALUE(symmetryY, "0")
        ATTRIBUTE_INSERT(symmetryY, "symmetryY:Dimension3")

        PROPERTY_DOUBLE(symmetryZ, "the direction of the positive symmetry axis, z component")
        ATTRIBUTE_QUANTITY(symmetryZ, "length")
        ATTRIBUTE_DEFAULT_VALUE(symmetryZ, "1")

        PROPERTY_DOUBLE(velocityX, "the bulk velocity of the point source, x component")
        ATTRIBUTE_QUANTITY(velocityX, "velocity")
        ATTRIBUTE_MIN_VALUE(velocityX, "[-100000 km/s")
        ATTRIBUTE_MAX_VALUE(velocityX, "100000 km/s]")
        ATTRIBUTE_DEFAULT_VALUE(velocityX, "0")
        ATTRIBUTE_RELEVANT_IF(velocityX, "Panchromatic")
        ATTRIBUTE_INSERT(velocityX, "Panchromatic&velocityX:Dimension3")

        PROPERTY_DOUBLE(velocityY, "the bulk velocity of the point source, y component")
        ATTRIBUTE_QUANTITY(velocityY, "velocity")
        ATTRIBUTE_MIN_VALUE(velocityY, "[-100000 km/s")
        ATTRIBUTE_MAX_VALUE(velocityY, "100000 km/s]")
        ATTRIBUTE_DEFAULT_VALUE(velocityY, "0")
        ATTRIBUTE_RELEVANT_IF(velocityY, "Panchromatic")
        ATTRIBUTE_INSERT(velocityY, "Panchromatic&velocityY:Dimension3")

        PROPERTY_DOUBLE(velocityZ, "the bulk velocity of the point source, z component")
        ATTRIBUTE_QUANTITY(velocityZ, "velocity")
        ATTRIBUTE_MIN_VALUE(velocityZ, "[-100000 km/s")
        ATTRIBUTE_MAX_VALUE(velocityZ, "100000 km/s]")
        ATTRIBUTE_DEFAULT_VALUE(velocityZ, "0")
        ATTRIBUTE_RELEVANT_IF(velocityZ, "Panchromatic")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function reads the user-provided input file and caches the relevant information. */
    void setupSelfBefore() override;

    /** This function warns the user if this source's intrinsic wavelength range does not fully
        cover the configured wavelength range. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the point source, which depends on the (lack of)
        symmetry of its geometry. */
    int dimension() const override;

    /** This function returns true if the velocity of the source is nonzero. */
    bool hasVelocity() const override;

    /** This function implements the VelocityInterface interface. It returns the bulk velocity
         of this source, as configured by the user. */
    Vec velocity() const override;

    /** This function returns the wavelength range for this source. Outside this range, all
        luminosities are zero. This source's wavelength range is determined as the intersection of
        the simulation's source wavelength range (obtained from the simulation configuration) and
        the intrinsic wavelength range of the spectrum read from the input file.

        This function implements the SourceWavelengthRangeInterface interface. */
    Range wavelengthRange() const override;

    /** This function returns the luminosity \f$L\f$ (i.e. radiative power) of the source
        integrated over the wavelength range of primary sources (configured for the source system
        as a whole). */
    double luminosity() const override;

    /** This function returns the specific luminosity \f$L_\lambda\f$ (i.e. radiative power per
        unit of wavelength) of the source at the specified wavelength, or zero if the wavelength is
        outside the wavelength range of primary sources (configured for the source system as a
        whole) or if the source simply does not emit at the wavelength. */
    double specificLuminosity(double wavelength) const override;

    /** This function causes the photon packet \em pp to be launched from the source using the
        given history index and luminosity contribution. The function samples a wavelength and a
        propagation direction from the spectral and angular dependencies derived from the input
        file, constructs proper objects for angular distribution, polarization profile and bulk
        velocity, and finally emits the photon packet with these properties. */
    void launch(PhotonPacket* pp, size_t historyIndex, double L) const override;

    //======================== Data Members ========================

private:
    // spectral data initialized during setup
    ContSED* _sed{nullptr};  // the SED averaged over all angles as derived from the input file

    // sampling parameters initialized during setup
    double _xi{0.};                                      // the wavelength bias fraction
    WavelengthDistribution* _biasDistribution{nullptr};  // the wavelength bias distribution

    // velocity data initialized during setup
    VelocityInterface* _bvi{nullptr};  // pointer to an object offering the velocity interface, if needed
};

//////////////////////////////////////////////////////////////////////

#endif
