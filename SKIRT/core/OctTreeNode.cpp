/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "OctTreeNode.hpp"

////////////////////////////////////////////////////////////////////

// macros for easily accessing a particular child
#define CHILD_0 childAt(0)
#define CHILD_1 childAt(1)
#define CHILD_2 childAt(2)
#define CHILD_3 childAt(3)
#define CHILD_4 childAt(4)
#define CHILD_5 childAt(5)
#define CHILD_6 childAt(6)
#define CHILD_7 childAt(7)

////////////////////////////////////////////////////////////////////

void OctTreeNode::createChildren(int id)
{
    Vec rc = center();
    addChild(new OctTreeNode(this, id++, Box(xmin(), ymin(), zmin(), rc.x(), rc.y(), rc.z())));
    addChild(new OctTreeNode(this, id++, Box(rc.x(), ymin(), zmin(), xmax(), rc.y(), rc.z())));
    addChild(new OctTreeNode(this, id++, Box(xmin(), rc.y(), zmin(), rc.x(), ymax(), rc.z())));
    addChild(new OctTreeNode(this, id++, Box(rc.x(), rc.y(), zmin(), xmax(), ymax(), rc.z())));
    addChild(new OctTreeNode(this, id++, Box(xmin(), ymin(), rc.z(), rc.x(), rc.y(), zmax())));
    addChild(new OctTreeNode(this, id++, Box(rc.x(), ymin(), rc.z(), xmax(), rc.y(), zmax())));
    addChild(new OctTreeNode(this, id++, Box(xmin(), rc.y(), rc.z(), rc.x(), ymax(), zmax())));
    addChild(new OctTreeNode(this, id++, Box(rc.x(), rc.y(), rc.z(), xmax(), ymax(), zmax())));
}

////////////////////////////////////////////////////////////////////

TreeNode* OctTreeNode::child(Vec r)
{
    Vec rc = CHILD_0->rmax();
    int l = (r.x() < rc.x() ? 0 : 1) + (r.y() < rc.y() ? 0 : 2) + (r.z() < rc.z() ? 0 : 4);
    return children()[l];
}

////////////////////////////////////////////////////////////////////

void OctTreeNode::addNeighbors()
{
    // if we don't have any children, we can't add neighbors
    if (isChildless()) return;

    // Internal neighbors: each of the 8 new children has 3 neighbors among its siblings
    makeNeighbors(FRONT, CHILD_0, CHILD_1);
    makeNeighbors(RIGHT, CHILD_0, CHILD_2);
    makeNeighbors(TOP, CHILD_0, CHILD_4);
    makeNeighbors(RIGHT, CHILD_1, CHILD_3);
    makeNeighbors(TOP, CHILD_1, CHILD_5);
    makeNeighbors(FRONT, CHILD_2, CHILD_3);
    makeNeighbors(TOP, CHILD_2, CHILD_6);
    makeNeighbors(TOP, CHILD_3, CHILD_7);
    makeNeighbors(FRONT, CHILD_4, CHILD_5);
    makeNeighbors(RIGHT, CHILD_4, CHILD_6);
    makeNeighbors(RIGHT, CHILD_5, CHILD_7);
    makeNeighbors(FRONT, CHILD_6, CHILD_7);

    // The point where this node is split into its children
    double xc = CHILD_0->xmax();
    double yc = CHILD_0->ymax();
    double zc = CHILD_0->zmax();

    // The BACK neighbors of this node
    {
        for (auto neighbor : neighbors(BACK))
        {
            neighbor->deleteNeighbor(FRONT, this);
            if (neighbor->ymin() <= yc && neighbor->zmin() <= zc) makeNeighbors(FRONT, neighbor, CHILD_0);
            if (neighbor->ymax() >= yc && neighbor->zmin() <= zc) makeNeighbors(FRONT, neighbor, CHILD_2);
            if (neighbor->ymin() <= yc && neighbor->zmax() >= zc) makeNeighbors(FRONT, neighbor, CHILD_4);
            if (neighbor->ymax() >= yc && neighbor->zmax() >= zc) makeNeighbors(FRONT, neighbor, CHILD_6);
        }
    }
    // The FRONT neighbors of this node
    {
        for (auto neighbor : neighbors(FRONT))
        {
            neighbor->deleteNeighbor(BACK, this);
            if (neighbor->ymin() <= yc && neighbor->zmin() <= zc) makeNeighbors(BACK, neighbor, CHILD_1);
            if (neighbor->ymax() >= yc && neighbor->zmin() <= zc) makeNeighbors(BACK, neighbor, CHILD_3);
            if (neighbor->ymin() <= yc && neighbor->zmax() >= zc) makeNeighbors(BACK, neighbor, CHILD_5);
            if (neighbor->ymax() >= yc && neighbor->zmax() >= zc) makeNeighbors(BACK, neighbor, CHILD_7);
        }
    }
    // The LEFT neighbors of this node
    {
        for (auto neighbor : neighbors(LEFT))
        {
            neighbor->deleteNeighbor(RIGHT, this);
            if (neighbor->xmin() <= xc && neighbor->zmin() <= zc) makeNeighbors(RIGHT, neighbor, CHILD_0);
            if (neighbor->xmax() >= xc && neighbor->zmin() <= zc) makeNeighbors(RIGHT, neighbor, CHILD_1);
            if (neighbor->xmin() <= xc && neighbor->zmax() >= zc) makeNeighbors(RIGHT, neighbor, CHILD_4);
            if (neighbor->xmax() >= xc && neighbor->zmax() >= zc) makeNeighbors(RIGHT, neighbor, CHILD_5);
        }
    }
    // The RIGHT neighbors of this node
    {
        for (auto neighbor : neighbors(RIGHT))
        {
            neighbor->deleteNeighbor(LEFT, this);
            if (neighbor->xmin() <= xc && neighbor->zmin() <= zc) makeNeighbors(LEFT, neighbor, CHILD_2);
            if (neighbor->xmax() >= xc && neighbor->zmin() <= zc) makeNeighbors(LEFT, neighbor, CHILD_3);
            if (neighbor->xmin() <= xc && neighbor->zmax() >= zc) makeNeighbors(LEFT, neighbor, CHILD_6);
            if (neighbor->xmax() >= xc && neighbor->zmax() >= zc) makeNeighbors(LEFT, neighbor, CHILD_7);
        }
    }
    // The BOTTOM neighbors of this node
    {
        for (auto neighbor : neighbors(BOTTOM))
        {
            neighbor->deleteNeighbor(TOP, this);
            if (neighbor->xmin() <= xc && neighbor->ymin() <= yc) makeNeighbors(TOP, neighbor, CHILD_0);
            if (neighbor->xmax() >= xc && neighbor->ymin() <= yc) makeNeighbors(TOP, neighbor, CHILD_1);
            if (neighbor->xmin() <= xc && neighbor->ymax() >= yc) makeNeighbors(TOP, neighbor, CHILD_2);
            if (neighbor->xmax() >= xc && neighbor->ymax() >= yc) makeNeighbors(TOP, neighbor, CHILD_3);
        }
    }
    // The TOP neighbors of this node
    {
        for (auto neighbor : neighbors(TOP))
        {
            neighbor->deleteNeighbor(BOTTOM, this);
            if (neighbor->xmin() <= xc && neighbor->ymin() <= yc) makeNeighbors(BOTTOM, neighbor, CHILD_4);
            if (neighbor->xmax() >= xc && neighbor->ymin() <= yc) makeNeighbors(BOTTOM, neighbor, CHILD_5);
            if (neighbor->xmin() <= xc && neighbor->ymax() >= yc) makeNeighbors(BOTTOM, neighbor, CHILD_6);
            if (neighbor->xmax() >= xc && neighbor->ymax() >= yc) makeNeighbors(BOTTOM, neighbor, CHILD_7);
        }
    }
}

////////////////////////////////////////////////////////////////////
