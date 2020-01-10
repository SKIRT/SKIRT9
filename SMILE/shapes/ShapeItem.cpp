/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ShapeItem.hpp"

////////////////////////////////////////////////////////////////////

void ShapeItem::paint()
{
    pushState();
    paintSelf();
    for (auto child : children())
    {
        auto c = dynamic_cast<ShapeItem*>(child);
        if (c) c->paint();
    }
    popState();
}

////////////////////////////////////////////////////////////////////

void ShapeItem::pushState()
{
    auto p = dynamic_cast<ShapeItem*>(parent());
    if (p) p->pushState();
}

////////////////////////////////////////////////////////////////////

void ShapeItem::popState()
{
    auto p = dynamic_cast<ShapeItem*>(parent());
    if (p) p->popState();
}

////////////////////////////////////////////////////////////////////

void ShapeItem::paintSelf() {}

////////////////////////////////////////////////////////////////////

void ShapeItem::setColor(double r, double g, double b)
{
    auto p = dynamic_cast<ShapeItem*>(parent());
    if (p) p->setColor(r, g, b);
}

////////////////////////////////////////////////////////////////////

void ShapeItem::setWidth(double w)
{
    auto p = dynamic_cast<ShapeItem*>(parent());
    if (p) p->setWidth(w);
}

////////////////////////////////////////////////////////////////////

void ShapeItem::drawLine(double x1, double y1, double x2, double y2)
{
    auto p = dynamic_cast<ShapeItem*>(parent());
    if (p) p->drawLine(x1, y1, x2, y2);
}

////////////////////////////////////////////////////////////////////
