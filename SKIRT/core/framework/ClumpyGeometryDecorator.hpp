/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CLUMPYGEOMETRYDECORATOR_HPP
#define CLUMPYGEOMETRYDECORATOR_HPP

#include "GenGeometry.hpp"
#include "SmoothingKernel.hpp"

////////////////////////////////////////////////////////////////////

/** The ClumpyGeometryDecorator class is a geometry decorator that adds clumpiness to any geometry.
    It basically assigns a fraction \f$f\f$ of the mass of the original geometry to compact clumps,
    which are distributed statistically according to the same distribution. The properties of a
    ClumpyGeometryDecorator object are a reference to the original Geometry object being decorated,
    and the characteristics that describe the clumpiness, i.e. the fraction \f$f\f$ of the mass
    locked in clumps, the total number \f$N\f$ of clumps, the scale radius \f$h\f$ of a single
    clump, and the kernel \f$W({\bf{r}},h)\f$ that describes the mass distribution of a single
    clump. If the original geometry is characterized by the density
    \f$\rho_{\text{orig}}({\bf{r}})\f$, the new, clumpy stellar geometry is described by \f[
    \rho({\bf{r}}) = (1-f)\, \rho_{\text{orig}} ({\bf{r}}) + \frac{f}{N} \sum_{i=1}^N
    W({\bf{r}}-{\bf{r}}_i,h). \f] where \f${\bf{r}}_i\f$ is the location of the centre of the
    \f$i\f$'th clump, each of them randomly drawn from the three-dimensional probability density
    \f$p({\bf{r}})\, {\text{d}}{\bf{r}} = \rho_{\text{orig}}({\bf{r}})\, {\text{d}}{\bf{r}}\f$.

    By default, this class uses the standard random generator also used by other classes during
    setup. Consecutive executions of the same ski file will produce the same clump positions (even
    if the simulation is configured to have multiple parallel execution threads or processes). On
    the other hand, multiple occurrences of the ClumpyGeometryDecorator in a given ski file will
    always produce a different set of clump positions, because consecutive portions of the
    pseudo-random sequence are being employed. While this is usually just fine, in some models one
    might want to line up, for example, the clumps in a medium distribution with those in a source
    distribution. Therefore, if a nonzero value is specified for the \em seed property, the clump
    positions are generated using a temporary random number generator initialized with that seed.
    Configuring the same seed for two ClumpyGeometryDecorator instances will line up the respective
    clump positions, assuming the underlying geometry and the number of clumps are identical. */
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

        PROPERTY_ITEM(smoothingKernel, SmoothingKernel,
                      "the smoothing kernel that describes the density of a single clump")
        ATTRIBUTE_DEFAULT_VALUE(smoothingKernel, "CubicSplineSmoothingKernel")
        ATTRIBUTE_DISPLAYED_IF(smoothingKernel, "Level2")

        PROPERTY_INT(seed, "the seed for the random clump position generator, or zero to use the default generator")
        ATTRIBUTE_MIN_VALUE(seed, "0")
        ATTRIBUTE_MAX_VALUE(seed, "1000000")
        ATTRIBUTE_DEFAULT_VALUE(seed, "0")
        ATTRIBUTE_DISPLAYED_IF(seed, "Level3")

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
