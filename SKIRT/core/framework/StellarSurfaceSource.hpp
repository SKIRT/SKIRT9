/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STELLARSURFACESOURCE_HPP
#define STELLARSURFACESOURCE_HPP

#include "CenteredSource.hpp"

////////////////////////////////////////////////////////////////////

/** The StellarSurfaceSource class represents the surface of a sphere of radius \f$r_*\f$ from
    which radiation escapes in the outward direction. The emissivity of this source is anisotropic.
    There is no radiation inwards and the intensity is proportional to \f$\cos\theta'\f$ in the
    outward hemisphere, where \f$\theta'\f$ is the angle between the direction and the normal on
    the surface. */
class StellarSurfaceSource : public CenteredSource
{
    ITEM_CONCRETE(StellarSurfaceSource, CenteredSource,
                  "a stellar surface source with an anisotropic outward radiation field")

        PROPERTY_DOUBLE(stellarRadius, "the stellar radius")
        ATTRIBUTE_QUANTITY(stellarRadius, "length")
        ATTRIBUTE_MIN_VALUE(stellarRadius, "]0")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the intrinsic dimension of the geometry, which is 1 in this case: the
        geometry is spherically symmetric, even with the anisotropic emissivity. */
    int intrinsicDimension() const override;

    /** This function causes the photon packet \em pp to be launched from the source using the
        given history index, wavelength, weighted luminosity contribution, and redshift interface
        (corresponding to the bulk velocity of the source).

        Firstly, the function generates a random outward launch position \f${\bf{r}}\f$ on the
        stellar surface by drawing a random direction \f${\bf{u}}\f$ on the unit sphere and
        transforming it to take into account the stellar radius and center position. The direction
        \f${\bf{u}}\f$ actually corresponds to the outward normal to the stellar surface at the
        position \f${\bf{r}}\f$.

        Secondly, the function generates a random direction \f${\bf{k}}\f$ appropriate for the
        anisotropic outward radiation field at the location \f${\bf{r}}\f$, or equivalently, for
        the normal \f${\bf{u}}\f$. In a spherical coordinate system with the Z' axis along the
        normal \f${\bf{u}}\f$, the probability distribution is given by \f[ p(\theta',\varphi')\,
        {\text{d}}\theta'\, {\text{d}}\varphi' = \frac{1}{\pi}\cos\theta'
        \sin\theta'\,{\text{d}}\theta'\, {\text{d}}\varphi' \qquad 0 \leq \theta' < \frac{\pi}{2}.
        \f] Random angles \f$\theta'\f$ and \f$\varphi'\f$ can be determined by taking two uniform
        deviates \f${\cal{X}}_1\f$ and \f${\cal{X}}_2\f$ and solving the two equations \f[
        \begin{split} {\cal{X}}_1 &= \int_0^{\theta'}
        2\sin\theta^{\prime\prime}\cos\theta^{\prime\prime}\, {\text{d}}\theta^{\prime\prime} \\
        {\cal{X}}_2 &= \int_0^{\varphi'} \frac{{\text{d}}\varphi^{\prime\prime}}{2\pi} \end{split}
        \f] for \f$\theta'\f$ and \f$\varphi'\f$. The solution is readily found, \f[ \begin{split}
        \theta' &= \arcsin \sqrt{{\cal{X}}_1} \\ \varphi' &= 2\pi\,{\cal{X}}_2. \end{split} \f]
        Once these values have been determined, we need to determine the coordinates of the vector
        \f${\bf{k}}\f$ in the standard coordinate system. The link between the standard XYZ
        coordinate system and the one with the Z' axis determined by the normal \f${\bf{u}}\f$ is a
        joint rotation over an angle \f$\varphi\f$ along the Z-axis, followed by a rotation over an
        angle \f$\theta\f$ over the Y'-axis, where \f$\theta\f$ and \f$\varphi\f$ are the polar
        angle and azimuth of the position vector \f${\bf{r}}\f$. The connection is hence \f[
        \begin{bmatrix} k'_x \\ k'_y \\ k'_z \end{bmatrix} = \begin{bmatrix} \cos\theta & 0 &
        -\sin\theta \\ 0 & 1 & 0 \\ \sin\theta & 0 & \cos\theta \end{bmatrix} \begin{bmatrix}
        \cos\phi & \sin\phi & 0 \\ -\sin\phi & \cos\phi & 0 \\ 0 & 0 & 1 \end{bmatrix}
        \begin{bmatrix} k_x \\ k_y \\ k_z \end{bmatrix} \f] or in the other direction \f[
        \begin{bmatrix} k_x \\ k_y \\ k_z \end{bmatrix} = \begin{bmatrix} \cos\phi & -\sin\phi & 0
        \\ \sin\phi & \cos\phi & 0 \\ 0 & 0 & 1 \end{bmatrix} \begin{bmatrix} \cos\theta & 0 &
        \sin\theta \\ 0 & 1 & 0 \\ -\sin\theta & 0 & \cos\theta \end{bmatrix} \begin{bmatrix} k'_x
        \\ k'_y \\ k'_z \end{bmatrix} \f]

        Thirdly, the function sets up a thread-local object that serves the angular distribution
        interface for a particular launch position. Configured with the appropriate launch position
        \f${\bf{u}}\f$, this object is handed to the photon packet being launched so that peel-off
        photon packets can retrieve the appropriate bias factor for the instrument direction. This
        works even if there are multiple sources of this type because each thread handles a single
        photon packet at a time.

        The local angular distribution object implements a function that returns the normalized
        probability for a direction \f${\bf{k}}\f$, given that the point of emission is defined by
        the normal \f${\bf{u}}\f$. The appropriate probability distribution is given by \f[
        p({\bf{k}})\, {\text{d}}\Omega = \begin{cases} 4\cos\theta'\, {\text{d}}\Omega & 0 \leq
        \theta' < \frac{\pi}{2} \\ 0 & \frac{\pi}{2} \leq \theta' < \pi \end{cases} \f] Here
        \f$\theta'\f$ is the angle between the direction \f${\bf{k}}\f$ and the outward normal
        \f${\bf{u}}\f$ on the stellar surface, i.e. \f[ \cos\theta' = {\bf{k} \cdot \bf{r}}. \f]
        This distribution is correctly normalized in the sense that \f[ \frac{1}{4\pi} \iint
        p({\bf{k}})\, {\text{d}}\Omega = 1. \f]

        Finally, the function launches the photon packet, passing it all of the above information.
        */
    void launchSpecialty(PhotonPacket* pp, size_t historyIndex, double lambda, double Lw,
                         VelocityInterface* bvi) const override;
};

////////////////////////////////////////////////////////////////////

#endif
