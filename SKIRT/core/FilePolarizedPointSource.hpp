/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILEPOLARIZEDPOINTSOURCE_HPP
#define FILEPOLARIZEDPOINTSOURCE_HPP

#include "Direction.hpp"
#include "LuminosityNormalization.hpp"
#include "Source.hpp"
#include "StoredTable.hpp"
#include "VelocityInterface.hpp"
class ContSED;

//////////////////////////////////////////////////////////////////////

/** FilePolarizedPointSource represents a primary source limited to a single point in space and
    emitting polarized radiation with an axisymmetric angular dependence. The Stokes vector
    components of the emitted radiation, as a function of wavelength and inclination angle cosine,
    are loaded from a user-provided file. The bolometric output is characterized by a
    LuminosityNormalization object configured in the ski file. The class furthermore offers
    user-configured properties to specify the position of the point source, the orientation of the
    symmetry axis of the angular dependence, and a "bulk" velocity.

    <em>Luminosity normalization</em>

    The luminosity normalization property of this class can be configured in the ski file as usual,
    keeping in mind that the given value refers to the (specific or wavelength-integrated)
    luminosity integrated over the unit sphere.

    <em>Inclination angle</em>

    For the purposes of this class, the inclination angle \f$\theta\f$ is defined as the angle
    between the user-configured symmetry axis \f$\bf{s}\f$ and the propagation direction
    \f$\bf{k}\f$ of the emitted photon packet. The inclination cosine values given in the input
    file must conform to this convention.

    In formula form, assuming that both directions \f$\bf{s}\f$ and \f$\bf{k}\f$ are given as unit
    vectors, the inclination angle cosine simply is \f$\cos\theta=\bf{s}\cdot\bf{k}\f$.

    <em>Reference direction</em>

    The components of a Stokes vector (and hence the polarization angle derived from it) are
    specified relative to some reference direction \f$\bf{d}\f$ which must be perpendicular to the
    propagation direction \f$\bf{k}\f$ of the emitted photon packet. For the purposes of this
    class, the reference direction is taken to be the orthogonal projection of the user-configured
    symmetry axis \f$\bf{s}\f$ on the plane perpendicular to the propagation direction
    \f$\bf{k}\f$. The Stokes vector values given in the input file must conform to this convention.

    In formula form, assuming that the propagation direction \f$\bf{k}\f$ is given as a unit
    vector, the (unnormalized) orthogonal projection of \f$\bf{s}\f$ on the plane with normal
    \f$\bf{k}\f$ is given by \f$\bf{d} = \bf{s} - (\bf{k}\cdot\bf{s})\bf{k}\f$. Rather than the
    reference direction itself, the implementation instead uses the normal \f$\bf{n}\f$ to the
    reference plane, which is defined as the plane through both the reference direction and the
    propagation direction. The normalized normal vector \f$\bf{n}\f$ is obtained from

    \f[ \bf{n} = \frac{\bf{d} \times \bf{k}}{||\bf{d} \times \bf{k}||} = \frac{\bf{s} \times
    \bf{k}}{||\bf{s} \times \bf{k}||} \f]

    after substitution of the above expression for \f$\bf{d}\f$. If the symmetry axis and the
    propogation direction are (anti)parallel, the normal vector \f$\bf{n}\f$ is undefined and the
    emitted radiation is taken to be unpolarized.

    <em>Input file format</em>

    As indicated in the introduction above, the user-provided input file tabulates the Stokes
    vector components \f$I,Q,U,V\f$ of the emitted radiation as a function of wavelength and
    inclination angle cosine. In other words, the input file represents a two-dimensional table
    with four values in each entry. All information is bundled in a single binary file using the
    SKIRT stored table format (see StoredTable).

    All quantities are specified in internal SKIRT units:

    - The inclination axis specifies cosine grid points in the range [-1, 1], in increasing order.
    This corresponds to angular values in the range [180 deg, 0 deg]; note the reverse order. In
    case the outer grid points are not equal to -1 and 1, the Stokes vector values for cosines
    outside of the grid are simply clamped to the corresponding outer grid point. This implies
    that, even if the emission pattern is symmetric around the 90 deg inclination, it is \em not
    sufficient to provide values for half of the inclination range. Instead, the reflected pattern
    must be included as well.

    - The wavelength axis specifies the wavelength grid points expressed in meter, in increasing
    order. For wavelengths outside of this grid, the emission is considered to be zero.

    - The Stokes vector components \f$I,Q,U,V\f$ are expressed in Watt per meter, with a scale
    factor that is abitrary but identical for the four components and across all Stokes vectors in
    the file. It is important to ensure that the quantities are expressed per unit of wavelength,
    as opposed to per unit of frequency or energy. Other than that, the precise units are
    irrelevant because of the arbitrary scaling.

    <em>Example python script</em>

    Below is an excerpt from a python script used to create an appropiate input file for this class
    from the information in a set of text files. In particular, the script shows how to create the
    file in SKIRT stored table format using the PTS function pts.storedtable.io.writeStoredTable.
    It also illustrates some of the conversions that may need to be performed.

\code{.py}
import numpy as np
from pts.storedtable.io import writeStoredTable

# read the inclination grid (in deg); range is 0..90 deg
theta = ...
numTheta = len(theta)

# read the energy grid (in keV)
E = ...
numE = len(E)

# read the Stokes vector components (in ph/s/keV) for each inclination bin
Ie = np.zeros((numTheta, numE))
Qe = np.zeros((numTheta, numE))
Ue = np.zeros((numTheta, numE))
for i in range(numTheta):
    Ie[i], Qe[i], Ue[i] = ...

# add symmetrical inclination range and reverse the order so that range becomes 180..0 deg
theta = np.concatenate((180 - theta, theta[::-1]))
Ie = np.concatenate((Ie, Ie[::-1]))
Qe = np.concatenate((Qe, Qe[::-1]))
Ue = np.concatenate((Ue, Ue[::-1]))

# convert inclination from deg to cosine in range -1..1
costheta = np.cos(theta * np.pi/180)

# convert energy in keV to wavelength in m
h = 6.62606957e-34
c = 2.99792458e8
q = 1.602176634e-19
lam = 1e-3 * h * c / q / E

# convert flux from per-energy flavor to per-wavelength flavor,
# with arbitrary normalization so that values are on the order of 1
Ilam = Ie / lam**3 / 1e35
Qlam = Qe / lam**3 / 1e35
Ulam = Ue / lam**3 / 1e35

# reverse spectral range so that wavelengths are in increasing order
lam = np.flip(lam)
Ilam = np.flip(Ilam, axis=1)
Qlam = np.flip(Qlam, axis=1)
Ulam = np.flip(Ulam, axis=1)

# create zero V component
Vlam = np.zeros_like(Ilam)

# write stored table
writeStoredTable(inputname + ".stab",
                 ['costheta','lambda'], ['1','m'], ['lin','log'], [costheta,lam],
                 ['I','Q','U','V'], ['W/m','W/m','W/m','W/m'], ['log','lin','lin','lin'], [Ilam,Qlam,Ulam,Vlam])
\endcode

    */
class FilePolarizedPointSource : public Source, public VelocityInterface
{
    ITEM_CONCRETE(FilePolarizedPointSource, Source, "a primary point source with a polarized spectrum read from file")
        ATTRIBUTE_TYPE_DISPLAYED_IF(FilePolarizedPointSource, "Level2")
        ATTRIBUTE_TYPE_INSERT(FilePolarizedPointSource, "Dimension2,ContSED")

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

    // structure combining the input tables representing the user-provided Stokes vector components
    struct Tables
    {
        StoredTable<2> I, Q, U, V;
    };

private:
    // user-provided Stokes vector components as a function of wavelength and inclination angle cosine
    Tables _tables;

    // spectral data initialized during setup
    ContSED* _sed{nullptr};  // the mean SED (averaged over the unit sphere) as derived from the input file

    // sampling parameters initialized during setup
    double _xi{0.};                                      // the wavelength bias fraction
    WavelengthDistribution* _biasDistribution{nullptr};  // the wavelength bias distribution

    // symmetry axis initialized during setup
    Direction _sym;  // direction of symmetry axis

    // velocity data initialized during setup
    VelocityInterface* _bvi{nullptr};  // pointer to an object offering the velocity interface, if needed
};

//////////////////////////////////////////////////////////////////////

#endif
