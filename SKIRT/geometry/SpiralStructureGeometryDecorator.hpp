/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPIRALSTRUCTUREGEOMETRYDECORATOR_HPP
#define SPIRALSTRUCTUREGEOMETRYDECORATOR_HPP

#include "AxGeometry.hpp"
#include "GenGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** The SpiralStructureGeometryDecorator class is a geometry decorator that adds spiral
    structure to any axisymmetric geometry. The spiral arm perturbation (with an
    arbitrary weight factor) is a logarithmic spiral arm pattern, based on the
    formulation of Schechtman-Rook et al. (2012, ApJ, 746, 70). The decorator
    basically alters the uniform distribution in azimuth (by definition, the
    density \f$\rho_{\text{ax}}\f$ of the original geometry is independent of
    \f$\phi\f$). In formula form, the density of the new geometry behaves as \f[
    \rho(R,\phi,z) = \rho_{\text{ax}}(R,z)\, \xi(R,\phi) \f]
    where \f$\xi(R,\phi)\f$ is a perturbation given by
    \f[ \xi(R,\phi) = (1-w) + w\, C_N \sin^{2N}
    \left[\frac{m}{2}\left( \frac{\ln (R/R_0)}{\tan p}-(\phi-\phi_0)\right) +
    \frac{\pi}{4} \right]. \f] Apart from the reference to the original
    geometry that is being decorated, the model contains six parameters: the number
    of spiral arms \f$m\f$, the pitch angle \f$p\f$, the spiral arm radius and
    phase zero-points \f$R_0\f$ and \f$\phi_0\f$, the spiral perturbation weight
    \f$w\f$, and the integer index \f$N>0\f$ that sets the arm-interarm size ratio
    (larger values of \f$N\f$ correspond to larger arm-interarm size ratios). The
    factor \f$C_N\f$ is not a free parameter, but a normalization factor that ensures
    that the total mass equals one, \f[ C_N = \frac{\sqrt{\pi}\,
    \Gamma(N+1)}{\Gamma(N+\tfrac12)}. \f] For \f$N=1\f$ the expression for the
    perturbation reduces to \f[ \xi(R,\phi) = 1 +
    w\, \sin \left[m \left( \frac{\ln (R/R_0)}{\tan p} - (\phi-\phi_0)\right)
    \right], \f] as in Misiriotis et al. (2000, A&A, 353, 117). Note that
    the parameters \f$R_0\f$ and \f$\phi_0\f$ in fact have the same effect (both of
    them add an offset to the spiral structure). In principle one of them could be
    suppressed, but it is confortable to include both of them. */
class SpiralStructureGeometryDecorator : public GenGeometry
{
    ITEM_CONCRETE(SpiralStructureGeometryDecorator, GenGeometry,
                  "a decorator that adds spiral structure to any axisymmetric geometry")

        PROPERTY_ITEM(geometry, AxGeometry, "the axisymmetric geometry to be decorated with spiral structure")

        PROPERTY_INT(numArms, "the number of spiral arms")
        ATTRIBUTE_MIN_VALUE(numArms, "1")
        ATTRIBUTE_MAX_VALUE(numArms, "100")
        ATTRIBUTE_DEFAULT_VALUE(numArms, "1")

        PROPERTY_DOUBLE(pitchAngle, "the pitch angle")
        ATTRIBUTE_QUANTITY(pitchAngle, "posangle")
        ATTRIBUTE_MIN_VALUE(pitchAngle, "]0 deg")
        ATTRIBUTE_MAX_VALUE(pitchAngle, "90 deg[")
        ATTRIBUTE_DEFAULT_VALUE(pitchAngle, "10 deg")

        PROPERTY_DOUBLE(radiusZeroPoint, "the radius zero-point")
        ATTRIBUTE_QUANTITY(radiusZeroPoint, "length")
        ATTRIBUTE_MIN_VALUE(radiusZeroPoint, "]0")

        PROPERTY_DOUBLE(phaseZeroPoint, "the phase zero-point")
        ATTRIBUTE_QUANTITY(phaseZeroPoint, "posangle")
        ATTRIBUTE_MIN_VALUE(phaseZeroPoint, "[0 deg")
        ATTRIBUTE_MAX_VALUE(phaseZeroPoint, "360 deg]")
        ATTRIBUTE_DEFAULT_VALUE(phaseZeroPoint, "0 deg")
        ATTRIBUTE_DISPLAYED_IF(phaseZeroPoint, "Level2")

        PROPERTY_DOUBLE(perturbationWeight, "the weight of the spiral perturbation")
        ATTRIBUTE_MIN_VALUE(perturbationWeight, "]0")
        ATTRIBUTE_MAX_VALUE(perturbationWeight, "1]")

        PROPERTY_INT(index, "the arm-interarm size ratio index")
        ATTRIBUTE_MIN_VALUE(index, "0")
        ATTRIBUTE_MAX_VALUE(index, "10")
        ATTRIBUTE_DEFAULT_VALUE(index, "1")
        ATTRIBUTE_DISPLAYED_IF(index, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function calculates some frequently used values. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho({\bf{r}})\f$ at the position \f${\bf{r}}\f$. It
        just implements the analytical formula. */
    double density(Position bfr) const override;

    /** This function generates a random position from the geometry by
        drawing a random point from the three-dimensional probability density
        \f$p({\bf{r}})\,{\text{d}}{\bf{r}} = \rho({\bf{r}})\, {\text{d}}{\bf{r}}\f$. We use a
        combination of the conditional distribution technique and the rejection technique. */
    Position generatePosition() const override;

    /** This function returns the surface mass density along the X-axis, i.e.
        the integration of the mass density along the entire X-axis, \f[
        \Sigma_X = \int_{-\infty}^\infty \rho(x,0,0)\, {\text{d}}x.\f] This integral cannot be
        calculated analytically, but when averaged over all lines-of-sight in the equatorial plane,
        the contribution of the spiral arm perturbation cancels out, and we recover the
        X-axis surface density of the corresponding unperturbed model. */
    double SigmaX() const override;

    /** This function returns the surface mass density along the Y-axis, i.e.
        the integration of the mass density along the entire Y-axis, \f[
        \Sigma_Y = \int_{-\infty}^\infty \rho(0,y,0)\, {\text{d}}y.\f] This integral cannot be
        calculated analytically, but when averaged over all lines-of-sight in the equatorial plane,
        the contribution of the spiral arm perturbation cancels out, and we recover the
        Y-axis surface density of the corresponding unperturbed model. */
    double SigmaY() const override;

    /** This function returns the surface mass density along the Z-axis, i.e.
        the integration of the mass density along the entire Z-axis, \f[
        \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\, {\text{d}}z.\f] For the present decorator, this integral
        is not really well defined, as the logarithmic spiral perturbation winds ever stronger when we get closer
        to the Z-axis. We use the Z-axis surface density of the corresponding unperturbed model. */
    double SigmaZ() const override;

private:
    /** This private function implements the analytical formula for the perturbation \f$\xi(R,\phi)\f$. */
    double perturbation(double R, double phi) const;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation and backwards compatibility
    const int& _m{_numArms};
    const double& _p{_pitchAngle};
    const double& _R0{_radiusZeroPoint};
    const double& _phi0{_phaseZeroPoint};
    const double& _w{_perturbationWeight};
    const int& _N{_index};

    // data members initialized during setup
    double _tanp{0.};
    double _cn{0.};
    double _c{0.};
};

////////////////////////////////////////////////////////////////////

#endif
