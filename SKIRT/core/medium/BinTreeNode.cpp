/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BinTreeNode.hpp"

////////////////////////////////////////////////////////////////////

// macros for easily accessing a particular child
#define CHILD_0 childAt(0)
#define CHILD_1 childAt(1)

////////////////////////////////////////////////////////////////////

// constants for splitting direction (not an enum because we use modulo arithmetic)
#define XDIR 0
#define YDIR 1
#define ZDIR 2

////////////////////////////////////////////////////////////////////

void BinTreeNode::createChildren(int id)
{
    switch (level() % 3)
    {
        case XDIR:
        {
            double xc = center().x();
            addChild(new BinTreeNode(this, id++, Box(xmin(), ymin(), zmin(), xc, ymax(), zmax())));
            addChild(new BinTreeNode(this, id++, Box(xc, ymin(), zmin(), xmax(), ymax(), zmax())));
        }
        break;
        case YDIR:
        {
            double yc = center().y();
            addChild(new BinTreeNode(this, id++, Box(xmin(), ymin(), zmin(), xmax(), yc, zmax())));
            addChild(new BinTreeNode(this, id++, Box(xmin(), yc, zmin(), xmax(), ymax(), zmax())));
        }
        break;
        case ZDIR:
        {
            double zc = center().z();
            addChild(new BinTreeNode(this, id++, Box(xmin(), ymin(), zmin(), xmax(), ymax(), zc)));
            addChild(new BinTreeNode(this, id++, Box(xmin(), ymin(), zc, xmax(), ymax(), zmax())));
        }
        break;
    }
}

////////////////////////////////////////////////////////////////////

TreeNode* BinTreeNode::child(Vec r)
{
    switch (level() % 3)
    {
        case XDIR: return r.x() < CHILD_0->xmax() ? CHILD_0 : CHILD_1;
        case YDIR: return r.y() < CHILD_0->ymax() ? CHILD_0 : CHILD_1;
        case ZDIR: return r.z() < CHILD_0->zmax() ? CHILD_0 : CHILD_1;
    }
    return nullptr;
}

////////////////////////////////////////////////////////////////////

void BinTreeNode::addNeighbors()
{
    // if we don't have any children, we can't add neighbors
    if (isChildless()) return;

    switch (level() % 3)
    {
        case XDIR:
        {
            double xc = CHILD_0->xmax();

            // Internal neighbors
            makeNeighbors(FRONT, CHILD_0, CHILD_1);

            // The BACK neighbors of this node
            {
                for (auto neighbor : neighbors(BACK))
                {
                    neighbor->deleteNeighbor(FRONT, this);
                    makeNeighbors(FRONT, neighbor, CHILD_0);
                }
            }
            // The FRONT neighbors of this node
            {
                for (auto neighbor : neighbors(FRONT))
                {
                    neighbor->deleteNeighbor(BACK, this);
                    makeNeighbors(BACK, neighbor, CHILD_1);
                }
            }
            // The LEFT neighbors of this node
            {
                for (auto neighbor : neighbors(LEFT))
                {
                    neighbor->deleteNeighbor(RIGHT, this);
                    if (neighbor->xmin() <= xc) makeNeighbors(RIGHT, neighbor, CHILD_0);
                    if (neighbor->xmax() >= xc) makeNeighbors(RIGHT, neighbor, CHILD_1);
                }
            }
            // The RIGHT neighbors of this node
            {
                for (auto neighbor : neighbors(RIGHT))
                {
                    neighbor->deleteNeighbor(LEFT, this);
                    if (neighbor->xmin() <= xc) makeNeighbors(LEFT, neighbor, CHILD_0);
                    if (neighbor->xmax() >= xc) makeNeighbors(LEFT, neighbor, CHILD_1);
                }
            }
            // The BOTTOM neighbors of this node
            {
                for (auto neighbor : neighbors(BOTTOM))
                {
                    neighbor->deleteNeighbor(TOP, this);
                    if (neighbor->xmin() <= xc) makeNeighbors(TOP, neighbor, CHILD_0);
                    if (neighbor->xmax() >= xc) makeNeighbors(TOP, neighbor, CHILD_1);
                }
            }
            // The TOP neighbors of this node
            {
                for (auto neighbor : neighbors(TOP))
                {
                    neighbor->deleteNeighbor(BOTTOM, this);
                    if (neighbor->xmin() <= xc) makeNeighbors(BOTTOM, neighbor, CHILD_0);
                    if (neighbor->xmax() >= xc) makeNeighbors(BOTTOM, neighbor, CHILD_1);
                }
            }
        }
        break;
        case YDIR:
        {
            double yc = CHILD_0->ymax();

            // Internal neighbors
            makeNeighbors(RIGHT, CHILD_0, CHILD_1);

            // The BACK neighbors of this node
            {
                for (auto neighbor : neighbors(BACK))
                {
                    neighbor->deleteNeighbor(FRONT, this);
                    if (neighbor->ymin() <= yc) makeNeighbors(FRONT, neighbor, CHILD_0);
                    if (neighbor->ymax() >= yc) makeNeighbors(FRONT, neighbor, CHILD_1);
                }
            }
            // The FRONT neighbors of this node
            {
                for (auto neighbor : neighbors(FRONT))
                {
                    neighbor->deleteNeighbor(BACK, this);
                    if (neighbor->ymin() <= yc) makeNeighbors(BACK, neighbor, CHILD_0);
                    if (neighbor->ymax() >= yc) makeNeighbors(BACK, neighbor, CHILD_1);
                }
            }
            // The LEFT neighbors of this node
            {
                for (auto neighbor : neighbors(LEFT))
                {
                    neighbor->deleteNeighbor(RIGHT, this);
                    makeNeighbors(RIGHT, neighbor, CHILD_0);
                }
            }
            // The RIGHT neighbors of this node
            {
                for (auto neighbor : neighbors(RIGHT))
                {
                    neighbor->deleteNeighbor(LEFT, this);
                    makeNeighbors(LEFT, neighbor, CHILD_1);
                }
            }
            // The BOTTOM neighbors of this node
            {
                for (auto neighbor : neighbors(BOTTOM))
                {
                    neighbor->deleteNeighbor(TOP, this);
                    if (neighbor->ymin() <= yc) makeNeighbors(TOP, neighbor, CHILD_0);
                    if (neighbor->ymax() >= yc) makeNeighbors(TOP, neighbor, CHILD_1);
                }
            }
            // The TOP neighbors of this node
            {
                for (auto neighbor : neighbors(TOP))
                {
                    neighbor->deleteNeighbor(BOTTOM, this);
                    if (neighbor->ymin() <= yc) makeNeighbors(BOTTOM, neighbor, CHILD_0);
                    if (neighbor->ymax() >= yc) makeNeighbors(BOTTOM, neighbor, CHILD_1);
                }
            }
        }
        break;
        case ZDIR:
        {
            double zc = CHILD_0->zmax();

            // Internal neighbors
            makeNeighbors(TOP, CHILD_0, CHILD_1);

            // The BACK neighbors of this node
            {
                for (auto neighbor : neighbors(BACK))
                {
                    neighbor->deleteNeighbor(BOTTOM, this);
                    neighbor->deleteNeighbor(FRONT, this);
                    if (neighbor->zmin() <= zc) makeNeighbors(FRONT, neighbor, CHILD_0);
                    if (neighbor->zmax() >= zc) makeNeighbors(FRONT, neighbor, CHILD_1);
                }
            }
            // The FRONT neighbors of this node
            {
                for (auto neighbor : neighbors(FRONT))
                {
                    neighbor->deleteNeighbor(BOTTOM, this);
                    neighbor->deleteNeighbor(BACK, this);
                    if (neighbor->zmin() <= zc) makeNeighbors(BACK, neighbor, CHILD_0);
                    if (neighbor->zmax() >= zc) makeNeighbors(BACK, neighbor, CHILD_1);
                }
            }
            // The LEFT neighbors of this node
            {
                for (auto neighbor : neighbors(LEFT))
                {
                    neighbor->deleteNeighbor(BOTTOM, this);
                    neighbor->deleteNeighbor(RIGHT, this);
                    if (neighbor->zmin() <= zc) makeNeighbors(RIGHT, neighbor, CHILD_0);
                    if (neighbor->zmax() >= zc) makeNeighbors(RIGHT, neighbor, CHILD_1);
                }
            }
            // The RIGHT neighbors of this node
            {
                for (auto neighbor : neighbors(RIGHT))
                {
                    neighbor->deleteNeighbor(BOTTOM, this);
                    neighbor->deleteNeighbor(LEFT, this);
                    if (neighbor->zmin() <= zc) makeNeighbors(LEFT, neighbor, CHILD_0);
                    if (neighbor->zmax() >= zc) makeNeighbors(LEFT, neighbor, CHILD_1);
                }
            }
            // The BOTTOM neighbors of this node
            {
                for (auto neighbor : neighbors(BOTTOM))
                {
                    neighbor->deleteNeighbor(BOTTOM, this);
                    neighbor->deleteNeighbor(TOP, this);
                    makeNeighbors(TOP, neighbor, CHILD_0);
                }
            }
            // The TOP neighbors of this node
            {
                for (auto neighbor : neighbors(TOP))
                {
                    neighbor->deleteNeighbor(BOTTOM, this);
                    neighbor->deleteNeighbor(BOTTOM, this);
                    makeNeighbors(BOTTOM, neighbor, CHILD_1);
                }
            }
        }
        break;
    }
}

////////////////////////////////////////////////////////////////////
