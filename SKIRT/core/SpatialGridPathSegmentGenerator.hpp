/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPATIALGRIDPATHSEGMENTGENERATOR_HPP
#define SPATIALGRIDPATHSEGMENTGENERATOR_HPP

#include "PathSegmentGenerator.hpp"
#include "SpatialGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** A SpatialGridPathSegmentGenerator object -- TO DO. */
class SpatialGridPathSegmentGenerator final
{
    // ------- Constructing and destructing -------

public:
    /** TO DO. */
    SpatialGridPathSegmentGenerator(const SpatialGridPath* path, const SpatialGrid* grid)
        : _generator(grid->createPathSegmentGenerator(path))
    {}

    /** TO DO. */
    ~SpatialGridPathSegmentGenerator() { delete _generator; }

    // ------- Generating and retrieving path segments -------

public:
    /** True if segment is available; false if no more segments. */
    bool next() { return _generator->next(); }

    /** TO DO. */
    int m() const { return _generator->m(); }

    /** TO DO. */
    double ds() const { return _generator->ds(); }

    // ------- Data members -------

private:
    PathSegmentGenerator* _generator{nullptr};
};

//////////////////////////////////////////////////////////////////////

#endif
