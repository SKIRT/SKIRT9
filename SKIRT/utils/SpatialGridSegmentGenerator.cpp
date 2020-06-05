/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory(), Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpatialGridSegmentGenerator.hpp"
#include "Box.hpp"

////////////////////////////////////////////////////////////////////

void SpatialGridSegmentGenerator::moveInside(const Box& box, double eps)
{
    // initialize to empty segment with zero length
    setEmptySegment();

    // --> x direction
    if (rx() <= box.xmin())
    {
        if (kx() <= 0.0)
        {
            setoutsidepos();
        }
        else
        {
            setEmptySegment((box.xmin() - rx()) / kx());
            setrx(box.xmin() + eps);
            propagatey();
            propagatez();
        }
    }
    else if (rx() >= box.xmax())
    {
        if (kx() >= 0.0)
        {
            setoutsidepos();
        }
        else
        {
            setEmptySegment((box.xmax() - rx()) / kx());
            setrx(box.xmax() - eps);
            propagatey();
            propagatez();
        }
    }

    // --> y direction
    if (ry() <= box.ymin())
    {
        if (ky() <= 0.0)
        {
            setoutsidepos();
        }
        else
        {
            setEmptySegment((box.ymin() - ry()) / ky());
            setry(box.ymin() + eps);
            propagatex();
            propagatez();
        }
    }
    else if (ry() >= box.ymax())
    {
        if (ky() >= 0.0)
        {
            setoutsidepos();
        }
        else
        {
            setEmptySegment((box.ymax() - ry()) / ky());
            setry(box.ymax() - eps);
            propagatex();
            propagatez();
        }
    }

    // --> z direction
    if (rz() <= box.zmin())
    {
        if (kz() <= 0.0)
        {
            setoutsidepos();
        }
        else
        {
            setEmptySegment((box.zmin() - rz()) / kz());
            setrz(box.zmin() + eps);
            propagatex();
            propagatey();
        }
    }
    else if (rz() >= box.zmax())
    {
        if (kz() >= 0.0)
        {
            setoutsidepos();
        }
        else
        {
            setEmptySegment((box.zmax() - rz()) / kz());
            setrz(box.zmax() - eps);
            propagatex();
            propagatey();
        }
    }
}

////////////////////////////////////////////////////////////////////
