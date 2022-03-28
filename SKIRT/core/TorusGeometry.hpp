/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TORUSGEOMETRY_HPP
#define TORUSGEOMETRY_HPP

#include "AxGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** The TorusGeometry class is a subclass of the AxGeometry class and describes the geometry of an
    axisymmetric torus as assumed to be present in the centre of active galactic nuclei (AGN). This
    geometry is described by a radial power-law density with a finite opening angle; see Stalevski
    et al. (2012, MNRAS, 420, 2756–2772) and Granato & Danese (1994, MNRAS, 268, 235). In
    formula, it is most easily expressed in spherical coordinates as \f[ \rho(r,\theta) = A\,
    r^{-p}\,{\text{e}}^{-q|\cos\theta|} \quad\text{for } r_{\text{min}}<r<r_{\text{max}} \text{ and
    } \frac{\pi}{2}-\Delta<\theta<\frac{\pi}{2} +\Delta. \f] There are five free parameters
    describing this dust geometry: the inner and outer radii \f$r_{\text{min}}\f$ and
    \f$r_{\text{max}}\f$ of the torus, the radial power law index \f$p\f$, the polar index \f$q\f$
    and the angle \f$\Delta\f$ describing the opening angle of the torus.

    If the dusty system under consideration is in the vicinity of an AGN central engine or another
    source which is luminous enough to heat the dust up to sublimation temperature, the inner
    radius should correspond to sublimation radius and scale as \f$ r_{\text{min}} \propto
    L(\theta)^{0.5}\f$ (Barvainis, 1987, ApJ, 320, 537, eq (5)). If the primary source assumes
    anisotropic emission, the inner radius must follow the same dependence as the distribution of
    the primary source luminosity. Otherwise, dust temperature on the inner boundary of geometry is
    very likely to be under- or over-estimated. Thus, if the NetzerAccretionDiskGeometry
    distribution is chosen to describe primary source emission, it is recommended to turn on the
    anisotropic inner radius option for the torus geometry. The inner radius will then be set by
    the following formula: \f[ r_{\text{min}} \propto (\cos\theta\,(2\cos\theta+1))^{0.5}.\f] This
    should allow dust to approach all the way to the primary central source in the equatorial
    plane. However, due to the finite resolution of dust cells, it may happen that some of the
    innermost cells end up with unphysically high temperatures. For this reason, there is an
    additional input parameter, the cutoff radius \f$r_{\text{cut}}\f$. The value of the cutoff
    radius is usually found after a few trial-and-error experiments by inspecting temperature
    distribution maps, until the inner wall of the geometry is at the expected sublimation
    temperature for a given dust population.

    The total dust mass of the model corresponds to the mass of the original geometry, before the
    inner wall is reshaped to account for anisotropy; the difference is usually rather small. */
class TorusGeometry : public AxGeometry
{
    ITEM_CONCRETE(TorusGeometry, AxGeometry, "a torus geometry")

        PROPERTY_DOUBLE(exponent, "the radial powerlaw exponent p of the torus")
        ATTRIBUTE_MIN_VALUE(exponent, "[0")

        PROPERTY_DOUBLE(index, "the polar index q of the torus")
        ATTRIBUTE_MIN_VALUE(index, "[0")

        PROPERTY_DOUBLE(openingAngle, "the half opening angle of the torus")
        ATTRIBUTE_QUANTITY(openingAngle, "posangle")
        ATTRIBUTE_MIN_VALUE(openingAngle, "[0 deg")
        ATTRIBUTE_MAX_VALUE(openingAngle, "90 deg]")

        PROPERTY_DOUBLE(minRadius, "the minimum radius of the torus")
        ATTRIBUTE_QUANTITY(minRadius, "length")
        ATTRIBUTE_MIN_VALUE(minRadius, "]0")

        PROPERTY_DOUBLE(maxRadius, "the maximum radius of the torus")
        ATTRIBUTE_QUANTITY(maxRadius, "length")
        ATTRIBUTE_MIN_VALUE(maxRadius, "]0")

        PROPERTY_BOOL(reshapeInnerRadius, "reshape the inner radius according to the Netzer luminosity profile")
        ATTRIBUTE_DEFAULT_VALUE(reshapeInnerRadius, "false")
        ATTRIBUTE_DISPLAYED_IF(reshapeInnerRadius, "Level2")

        PROPERTY_DOUBLE(cutoffRadius, "the inner cutoff radius of the torus")
        ATTRIBUTE_QUANTITY(cutoffRadius, "length")
        ATTRIBUTE_MIN_VALUE(cutoffRadius, "[0")
        ATTRIBUTE_DEFAULT_VALUE(cutoffRadius, "0")
        ATTRIBUTE_RELEVANT_IF(cutoffRadius, "reshapeInnerRadius")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function calculates some frequently used values. The normalization parameter \f$A\f$
        is set by the normalization condition that total mass equals one, i.e. \f[ 1 = 2\pi\, A\,
        \int_{\pi/2-\Delta}^{\pi/2+\Delta} e^{-q|\cos\theta|}\sin\theta\, {\text{d}}\theta
        \int_{r_{\text{min}}}^{r_{\text{max}}} r^{2-p}\, {\text{d}}r. \f] This results in \f[ A =
        \frac{q}{4\pi\, (1-{\text{e}}^{-q\sin\Delta})}\, \frac{1}{ {\text{gln}}_{p-2}\,
        r_{\text{max}} - {\text{gln}}_{p-2}\, r_{\text{min}} }, \f] with \f${\text{gln}}_p\, x\f$
        the generalized logarithm defined in SpecialFunctions::gln. If \f$q=0\f$, this expression
        reduces to \f[ A = \frac{1}{4\pi\,\sin\Delta\, ({\text{gln}}_{p-2}\, r_{\text{max}} -
        {\text{gln}}_{p-2}\, r_{\text{min}} )}. \f] */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho(R,z)\f$ at the cylindrical radius \f$R\f$ and
        height \f$z\f$. It just implements the analytical formula. */
    double density(double R, double z) const override;

    /** This function generates a random position from the torus geometry, by drawing a random
        point from the three-dimensional probability density \f$p({\bf{r}})\, {\text{d}}{\bf{r}} =
        \rho({\bf{r}})\, {\text{d}}{\bf{r}}\f$. For the torus geometry, the density is a separable
        function of \f$r\f$ and \f$\theta\f$, so that a random position can hence be constructed by
        combining random spherical coordinates, each chosen from their own probability
        distributions. A random azimuth \f$\phi\f$ is readily found by chosing a random deviate
        \f${\cal{X}}\f$ and setting \f$ \phi = 2\pi {\cal{X}} \f$.

        For the radial coordinate, the appropriate probability distribution is \f$
        p(r)\,{\text{d}}r \propto r^{2-p}\,{\text{d}}r \f$. A random radius is generated by picking
        a new uniform deviate \f${\cal{X}}\f$, and solving the equation \f[ {\cal{X}} =
        \int_{r_\text{min}}^r p(r')\, {\text{d}}r' \f] for \f$r\f$. For \f$p\ne3\f$ we find \f[
        {\cal{X}} = \frac{r^{3-p}-r_{\text{min}}^{3-p}}
        {r_{\text{max}}^{3-p}-r_{\text{min}}^{3-p}}. \f] Inverting this results in \f[ r = \left[
        (1-{\cal{X}})\,r_{\text{min}}^{3-p} + {\cal{X}}\,r_{\text{max}}^{3-p}
        \right]^{\frac{1}{3-p}}. \f] For \f$p=3\f$ this expression does not hold, and for
        \f$p\approx3\f$ it breaks down numerically. So for \f$p\approx3\f$ we can write the general
        expression \f[ r = {\text{gexp}}_{p-2} \Big[ {\text{gln}}_{p-2}\, r_{\text{min}} +
        {\cal{X}}\,( {\text{gln}}_{p-2}\, r_{\text{max}} - {\text{gln}}_{p-2}\, r_{\text{min}} )
        \Bigr]. \f] In this expression, \f${\text{gln}}_p\,x\f$ and \f${\text{gexp}}_p\,x\f$ are
        the generalized logarithm and exponential functions defined in SpecialFunctions::gln and
        SpecialFunctions::gexp respectively.

        Finally, for the polar angle, the appropriate distribution function is \f[ p(\theta)\,
        {\text{d}}\theta \propto e^{-q|\cos\theta|}\sin\theta\, {\text{d}}\theta. \f] A random
        polar angle is generated by picking a new uniform deviate \f${\cal{X}}\f$, and solving the
        equation \f[ {\cal{X}} = \int_0^\theta p(\theta')\, {\text{d}}\theta' \f] for \f$\theta\f$.
        We obtain after some calculation \f[ {\cal{X}} = \begin{cases} \; \dfrac12 \left( 1 -
        \dfrac{1-{\text{e}}^{-q\cos\theta}}{1-{\text{e}}^{-q\sin\Delta}} \right) & \quad\text{for }
        \frac{\pi}{2}-\Delta < \theta < \frac{\pi}{2} \\[1.2em] \;\dfrac12 \left( 1 +
        \dfrac{1-{\text{e}}^{q\cos\theta}}{1-{\text{e}}^{-q\sin\Delta}} \right) & \quad\text{for }
        \frac{\pi}{2} < \theta < \frac{\pi}{2}+\Delta \end{cases} \f] Inverting this gives \f[
        \cos\theta = \begin{cases}\; -\dfrac{1}{q} \ln\left[ 1-\left(1-
        {\text{e}}^{-q\sin\Delta}\right) (1-2{\cal{X}}) \right] & \quad\text{if
        $0<{\cal{X}}<\tfrac12$} \\[1.2em] \; \dfrac{1}{q} \ln\left[ 1-\left(1
        -{\text{e}}^{-q\sin\Delta}\right) (2{\cal{X}}-1) \right] & \quad\text{if
        $\tfrac12<{\cal{X}}<1$} \end{cases}. \f] */
    Position generatePosition() const override;

    /** This function returns the radial surface density, i.e. the integration of the density along
        a line in the equatorial plane starting at the centre of the coordinate system, \f[
        \Sigma_R = \int_0^\infty \rho(R,0)\,{\text{d}}R. \f] For the torus geometry, \f[ \Sigma_R =
        A\, ( {\text{gln}}_p\, r_{\text{max}} - {\text{gln}}_p\, r_{\text{min}} ) \f] with
        \f${\text{gln}}_p\,x\f$ the generalized logarithm defined in SpecialFunctions::gln. */
    double SigmaR() const override;

    /** This function returns the Z-axis surface density, i.e. the integration of the density along
        the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\, {\text{d}}z. \f] For
        the torus geometry this integral is simply zero (we exclude the special limiting case where
        \f$\Delta=\tfrac{\pi}{2}\f$). */
    double SigmaZ() const override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation and backwards compatibility
    const double& _p{_exponent};
    const double& _q{_index};
    const double& _Delta{_openingAngle};
    const double& _rmin{_minRadius};
    const double& _rmax{_maxRadius};
    const bool& _rani{_reshapeInnerRadius};
    const double& _rcut{_cutoffRadius};

    // data members initialized during setup
    double _sinDelta{0.};
    double _smin{0.};
    double _sdiff{0.};
    double _tmin{0.};
    double _tmax{0.};
    double _A{0.};
};

////////////////////////////////////////////////////////////////////

#endif
