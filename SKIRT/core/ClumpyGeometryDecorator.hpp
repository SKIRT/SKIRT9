/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CLUMPYGEOMETRYDECORATOR_HPP
#define CLUMPYGEOMETRYDECORATOR_HPP

#include "GenGeometry.hpp"
#include "SmoothingKernel.hpp"

////////////////////////////////////////////////////////////////////

/** The ClumpyGeometryDecorator class is a geometry decorator that adds clumpiness to
    any  geometry. It basically assigns a fraction \f$f\f$ of the mass of the original
    geometry to compact clumps, which are distributed statistically according to the same
    distribution. The properties of a ClumpyGeometryDecorator object are a reference
    to the original Geometry object being decorated, and the characteristics that
    describe the clumpiness, i.e. the fraction \f$f\f$ of the mass locked in clumps, the total
    number \f$N\f$ of clumps, the scale radius \f$h\f$ of a single clump, and the kernel
    \f$W({\bf{r}},h)\f$ that describes the mass distribution of a single clump. If
    the original geometry is characterized by the density \f$\rho_{\text{orig}}({\bf{r}})\f$, the
    new, clumpy stellar geometry is described by \f[ \rho({\bf{r}}) = (1-f)\, \rho_{\text{orig}}
    ({\bf{r}}) + \frac{f}{N} \sum_{i=1}^N W({\bf{r}}-{\bf{r}}_i,h). \f] where \f${\bf{r}}_i\f$ is
    the location of the centre of the \f$i\f$'th clump, each of them drawn stochastically from the
    three-dimensional probability density \f$p({\bf{r}})\, {\text{d}}{\bf{r}} =
    \rho_{\text{orig}}({\bf{r}})\, {\text{d}}{\bf{r}}\f$.*/
class ClumpyGeometryDecorator : public GenGeometry
{
    ITEM_CONCRETE(ClumpyGeometryDecorator, GenGeometry, "a decorator that adds clumpiness to any geometry")

    PROPERTY_ITEM(geometry, Geometry, "the geometry to be made clumpy")

    PROPERTY_DOUBLE(clumpFraction, "the fraction of the mass locked up in clumps")
        ATTRIBUTE_MIN_VALUE(clumpFraction, "[0")
        ATTRIBUTE_MAX_VALUE(clumpFraction, "1]")

    PROPERTY_INT(numClumps, "the total number of clumps")
        ATTRIBUTE_MIN_VALUE(numClumps, "1")

    PROPERTY_DOUBLE(clumpRadius, "the scale radius of a single clump")
        ATTRIBUTE_QUANTITY(clumpRadius, "length")
        ATTRIBUTE_MIN_VALUE(clumpRadius, "]0")

    PROPERTY_BOOL(cutoffClumps, "cut off clumps at the boundary of the underlying geometry")
        ATTRIBUTE_DEFAULT_VALUE(cutoffClumps, "false")
        ATTRIBUTE_DISPLAYED_IF(cutoffClumps, "Level2")

    PROPERTY_ITEM(smoothingKernel, SmoothingKernel, "the smoothing kernel that describes the density of a single clump")
        ATTRIBUTE_DEFAULT_VALUE(smoothingKernel, "CubicSplineSmoothingKernel")
        ATTRIBUTE_DISPLAYED_IF(smoothingKernel, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function generates the \f$N\f$ random positions corresponding
        to the centers of the individual clumps. They are chosen as random positions
        generated from the original geometry that is being decorated. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho({\bf{r}})\f$ at the position
        \f${\bf{r}}\f$. */
    double density(Position bfr) const override;

    /** This function generates a random position from the geometry, by drawing a random
        point from the three-dimensional probability density \f$p({\bf{r}})\, {\text{d}}{\bf{r}} =
        \rho({\bf{r}})\, {\text{d}}{\bf{r}}\f$. */
    Position generatePosition() const override;

    /** This pure virtual function returns the X-axis surface density. It simply passes on the
        value returned by the geometry being decorated. */
    double SigmaX() const override;

    /** This pure virtual function returns the Y-axis surface density. It simply passes on the
        value returned by the geometry being decorated. */
    double SigmaY() const override;

    /** This pure virtual function returns the Z-axis surface density. It simply passes on the
        value returned by the geometry being decorated. */
    double SigmaZ() const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    std::vector<Vec> _clumpv;
};

////////////////////////////////////////////////////////////////////

#endif
