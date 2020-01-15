/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CUBICALBACKGROUNDSOURCE_HPP
#define CUBICALBACKGROUNDSOURCE_HPP

#include "CenteredSource.hpp"

////////////////////////////////////////////////////////////////////

/** The CubicalBackgroundSource class represents the surface of a cube (i.e. the combination of the
    six walls) from which radiation escapes in the inward direction. The emissivity is anisotropic:
    there is no radiation outwards and the inward emissivity from each wall is proportional to
    \f$\cos\theta'\f$, where \f$\theta'\f$ is the angle between the direction and the normal on the
    wall. */
class CubicalBackgroundSource : public CenteredSource
{
    ITEM_CONCRETE(CubicalBackgroundSource, CenteredSource,
                  "a cubical background source with an anisotropic inward radiation field")
        ATTRIBUTE_TYPE_INSERT(CubicalBackgroundSource, "Dimension3")

        PROPERTY_DOUBLE(edgeLength, "the edge length of the background cube")
        ATTRIBUTE_QUANTITY(edgeLength, "length")
        ATTRIBUTE_MIN_VALUE(edgeLength, "]0")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the intrinsic dimension of the geometry, which is 3 in this case. */
    int intrinsicDimension() const override;

    /** This function causes the photon packet \em pp to be launched from the source using the
        given history index, wavelength, weighted luminosity contribution, and redshift interface
        (corresponding to the bulk velocity of the source).

        Firstly, the function generates a random inward launch position \f${\bf{r}}\f$ on the
        background surface. To accomplish this, it randomly selects one of the six cube walls
        (represented by its outward normal \f${\bf{u}}\f$), generates a random position within the
        selected wall, and transforms this position to take into account the background cube's size
        and center position.

        Secondly, the function generates a random direction \f${\bf{k}}\f$ appropriate for the
        anisotropic inward radiation field at the location \f${\bf{r}}\f$, or equivalently, for the
        wall normal \f${\bf{u}}\f$.

        In a spherical coordinate system with the Z' axis normal to the surface at the position
        \f${\bf{r}}\f$, the probability distribution is given by \f[ p(\theta',\varphi')\,
        {\text{d}}\theta'\, {\text{d}}\varphi' = -\frac{1}{\pi}\cos\theta'
        \sin\theta'\,{\text{d}}\theta'\, {\text{d}}\varphi' \qquad \frac{\pi}{2} \leq \theta' <
        \pi. \f] Random angles \f$\theta'\f$ and \f$\varphi'\f$ can be determined by taking two
        uniform deviates \f${\cal{X}}_1\f$ and \f${\cal{X}}_2\f$ and solving the two equations \f[
        \begin{split} {\cal{X}}_1 &= -\int_{\pi/2}^{\theta'}
        2\sin\theta^{\prime\prime}\cos\theta^{\prime\prime}\,{\text{d}}\theta^{\prime\prime} \\
        {\cal{X}}_2 &= \int_0^{\varphi'} \frac{{\text{d}}\varphi^{\prime\prime}}{2\pi} \end{split}
        \f] for \f$\theta'\f$ and \f$\varphi'\f$. The solution is readily found, \f[ \begin{split}
        \theta' &= \pi - \arccos \sqrt{{\cal{X}}_1} \\ \varphi' &= 2\pi\,{\cal{X}}_2. \end{split}
        \f] Once these values have been determined, we need to rotate such that the Z' axis is
        perpendicular to the wall corresponding to the position \f${\bf{r}}\f$.

        Thirdly, the function sets up a thread-local object that serves the angular distribution
        interface for a particular launch position. Configured with the appropriate launch wall
        normal \f${\bf{u}}\f$, this object is handed to the photon packet being launched so that
        peel-off photon packets can retrieve the appropriate bias factor for the instrument
        direction. This works even if there are multiple sources of this type because each thread
        handles a single photon packet at a time.

        The local angular distribution object implements a function that returns the normalized
        probability for a direction \f${\bf{k}}\f$, given that the point of emission is defined by
        the normal \f${\bf{u}}\f$. The appropriate probability distribution is given by \f[
        p({\bf{k}})\, {\text{d}}\Omega = \begin{cases} 0 & 0 \leq \theta' < \frac{\pi}{2} \\
        -4\cos\theta'\, {\text{d}}\Omega & \frac{\pi}{2} \leq \theta' < \pi \end{cases} \f] Here
        \f$\theta'\f$ is the angle between the direction \f${\bf{k}}\f$ and the outward normal
        \f${\bf{u}}\f$ of the launch wall, i.e. \f[ \cos\theta' = {\bf{k} \cdot \bf{r}}. \f] This
        distribution is correctly normalized in the sense that \f[ \frac{1}{4\pi} \iint
        p({\bf{k}})\, {\text{d}}\Omega = 1. \f]

        Finally, the function launches the photon packet, passing it all of the above information.
        */
    void launchSpecialty(PhotonPacket* pp, size_t historyIndex, double lambda, double Lw,
                         VelocityInterface* bvi) const override;
};

////////////////////////////////////////////////////////////////////

#endif
