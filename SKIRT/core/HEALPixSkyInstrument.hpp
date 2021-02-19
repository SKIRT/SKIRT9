/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef HEALPIXSKYINSTRUMENT_HPP
#define HEALPIXSKYINSTRUMENT_HPP

#include "AllSkyProjection.hpp"
#include "HomogeneousTransform.hpp"
#include "Instrument.hpp"

////////////////////////////////////////////////////////////////////

/** The HEALPixSkyInstrument class provides an all-sky instrument that does not perform any
    projection, but instead records fluxes onto a HEALPix tessellation of the sky sphere that
    guarantees equal surface areas for all pixels. (See the AllSkyInstrument class for an
    instrument that does perform a projection).

    This class is based on the official HEALPix algorithm by Górski, K. M. et al. (2005)
    (https://ui.adsabs.harvard.edu/abs/2005ApJ...622..759G/abstract), as available from
    https://healpix.sourceforge.io, but only uses the relevant part of the algorithm.

    In short, the HEALPix tessellation starts from a base tessellation of the sphere into 12
    quadrilaterals. 4 quadrilaterals are formed by subdividing the equator into 4 equal segments,
    so that each quadrilateral has two equatorial segment endpoints as vertices, complemented with
    two vertices that are obtained by translating the segment midpoint away from the equator over a
    zenith angle \f$\theta{}=\arccos(2/3)\f$ in both directions. The remaining 8 quadrilaterals are
    formed by connecting the north or south pole of the sphere with 3 vertices of already existing
    quadrilaterals.

    Within this base tessellation, higher order HEALPix tessellations are constructed by
    subdividing the 12 base quadrilaterals, as if these were simple squares in a flat Euclidian
    space. To guarantee that all pixels have the same surface area, this subdivision needs to use
    the same number of sub-pixels along each axis of the base quadrilateral. This number of
    sub-pixels is referred to as the HEALPix resolution parameter \f$N_\mathrm{side}\f$. It is
    common practice within the community to restrict \f$N_\mathrm{side}=2^k\f$, with \f$k\f$ an
    integer order, although this is not strictly required. We will use this parametrisation, as it
    allows for a more efficient pixel location algorithm.

    After subdivision, the HEALPix tessellation will contain \f$12N_\mathrm{side}^2\f$ pixels.
    These pixels are ordered so that their centres lie on \f$4N_\mathrm{side}-1\f$ rings of equal
    latitude. On each of these rings, pixel centres are equidistantly spaced in longitude, but the
    number of pixels in each ring, \f$i\f$, depends on the latitude. For pixels with a ring index
    \f$j\in{}[N_\mathrm{side},3N_\mathrm{side}]\f$, the number of pixels in the ring is constant
    and equal to \f$i=4N_\mathrm{side}\f$. Within the polar rings (\f$j < N_\mathrm{side}\f$ or
    \f$j>3N_\mathrm{side}\f$), the number of pixels decreases by 4 for every ring closer to the
    pole. For the polar rings in the northern hemisphere, the number of pixels is equal to
    \f$i=4j\f$. For the southern hemisphere, it is equal to \f$i=4(N_\mathrm{side}-j)\f$.

    HEALPix pixels can be ordered in two ways, depending on the application of interest. For many
    analysis tools, a hierarchical ordering is advantageous, whereby the index of neighbouring
    pixels are related. In practice, this so called \em nested ordering subdivides every base
    HEALPix pixel into a hierarchical quadtree structure. This is the approach used by e.g. the
    public Planck data maps. For visualisation purposes, the so called \em ring ordering is more
    appropriate, whereby the pixels are indexed per ring from northern to southern hemisphere, and
    with increasing longitude along each ring.

    For our implementation, we use an alternative ordering that closely matches the *ring*
    ordering. We compute the indices \f$j\f$ and \f$i\f$ for each pixel as if we were going to
    compute the *ring* pixel index, but then use these indices as the 2D pixel coordinates in a
    \f$(4N_\mathrm{side}-1)\times{}4N_\mathrm{side}\f$ image. Evidently, this leads to an excess of
    \f$4N_\mathrm{side}^2-16N\f$ pixels in the polar regions, similar to the empty pixels in
    projected AllSkyInstrument images. The advantage of this approach is that we can easily use the
    existing 2D image functionality for output, and that we can easily map output pixels to the
    corresponding HEALPix pixels in ring ordering during analysis.

    An HEALPix sky instrument consists of a sphere with a given radius, centered at the location of
    the instrument. Only photon packets originating outside of this sphere are detected. When a
    photon packet arrives, the spherical coordinates of its originating position relative to the
    instrument coordinate system are determined. The two resulting angular coordinates are
    transformed to pixel coordinates on the HEALPix sphere. The radial coordinate (i.e. the
    distance from the photon packet's origin to the instrument location) is used for adjusting the
    photon packet's luminosity contribution to the surface brightness in its target pixel. More
    details on the flux calibration for this "local" instrument are provided in the header of the
    FluxRecorder class.

    This instrument does \em not support recording of spatially integrated flux densities. */
class HEALPixSkyInstrument : public Instrument
{
    ITEM_CONCRETE(HEALPixSkyInstrument, Instrument,
                  "a HEALPix all-sky instrument (for Planck-like observations inside a model)")
        ATTRIBUTE_TYPE_DISPLAYED_IF(HEALPixSkyInstrument, "Level2")

        PROPERTY_INT(order, "HEALPix order")
        ATTRIBUTE_MIN_VALUE(order, "0")
        ATTRIBUTE_MAX_VALUE(order, "15")
        ATTRIBUTE_DEFAULT_VALUE(order, "6")

        PROPERTY_DOUBLE(radius, "the radius of the observer's all-sky sphere")
        ATTRIBUTE_QUANTITY(radius, "length")
        ATTRIBUTE_MIN_VALUE(radius, "]0")

        PROPERTY_DOUBLE(observerX, "the position of the observer, x component")
        ATTRIBUTE_QUANTITY(observerX, "length")

        PROPERTY_DOUBLE(observerY, "the position of the observer, y component")
        ATTRIBUTE_QUANTITY(observerY, "length")

        PROPERTY_DOUBLE(observerZ, "the position of the observer, z component")
        ATTRIBUTE_QUANTITY(observerZ, "length")

        PROPERTY_DOUBLE(crossX, "the position of the crosshair, x component")
        ATTRIBUTE_QUANTITY(crossX, "length")
        ATTRIBUTE_DEFAULT_VALUE(crossX, "0")

        PROPERTY_DOUBLE(crossY, "the position of the crosshair, y component")
        ATTRIBUTE_QUANTITY(crossY, "length")
        ATTRIBUTE_DEFAULT_VALUE(crossY, "0")

        PROPERTY_DOUBLE(crossZ, "the position of the crosshair, z component")
        ATTRIBUTE_QUANTITY(crossZ, "length")
        ATTRIBUTE_DEFAULT_VALUE(crossZ, "0")

        PROPERTY_DOUBLE(upX, "the upwards direction, x component")
        ATTRIBUTE_QUANTITY(upX, "length")
        ATTRIBUTE_DEFAULT_VALUE(upX, "0")

        PROPERTY_DOUBLE(upY, "the upwards direction, y component")
        ATTRIBUTE_QUANTITY(upY, "length")
        ATTRIBUTE_DEFAULT_VALUE(upY, "0")

        PROPERTY_DOUBLE(upZ, "the upwards direction, z component")
        ATTRIBUTE_QUANTITY(upZ, "length")
        ATTRIBUTE_DEFAULT_VALUE(upZ, "1")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that all attribute values have been appropriately set and performs
        setup for the instrument. Its most important task is to determine an appropriate
        transformation given the instrument position and configuration. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function determines whether the specified instrument has the same observer type,
        position and viewing direction as the receiving instrument, and if so, calls the
        setSameObserverAsPreceding() function to remember the fact. */
    void determineSameObserverAsPreceding(const Instrument* precedingInstrument) override;

    /** Returns the direction towards the observer from the given photon packet launching
        position, expressed in model coordinates. */
    Direction bfkobs(const Position& bfr) const override;

    /** Returns the direction along the positive y-axis in a plane normal to the vector towards the
        observer from the given photon packet launching position, expressed in model coordinates.
        */
    Direction bfky(const Position& bfr) const override;

protected:
    /** This function simulates the detection of a photon packet by the instrument. */
    void detect(PhotonPacket* pp) override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation
    const double& _Ox{_observerX};
    const double& _Oy{_observerY};
    const double& _Oz{_observerZ};
    const double& _Cx{_crossX};
    const double& _Cy{_crossY};
    const double& _Cz{_crossZ};
    const double& _Ux{_upX};
    const double& _Uy{_upY};
    const double& _Uz{_upZ};

    // data members derived from the published attributes during setup
    int _Nx{0};                       // number of pixels in the x direction
    int _Ny{0};                       // number of pixels in the y direction
    int _Nside{0};                    // number of pixels in one direction of a HEALPix base pixel
    HomogeneousTransform _transform;  // transform from world to observer coordinates
};

////////////////////////////////////////////////////////////////////

#endif
