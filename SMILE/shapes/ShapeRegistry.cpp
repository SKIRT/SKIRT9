/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ShapeRegistry.hpp"
#include "BuildInfo.hpp"
#include "ColorDecorator.hpp"
#include "ItemRegistry.hpp"
#include "PolygonShape.hpp"
#include "RectangleShape.hpp"
#include "ShapeCanvas.hpp"
#include "ShapeGroup.hpp"
#include "ShapesUnitDef.hpp"
#include "WidthDecorator.hpp"

////////////////////////////////////////////////////////////////////

ShapeRegistry::ShapeRegistry()
{
    string version = BuildInfo::projectVersion();
    ItemRegistry::beginSchema("Shapes", "a Shapes parameter file", version, "shapes", "shapes-definition",
                              "ShapeCanvas", version, "");

    ItemRegistry::add<ShapeItem>();
    ItemRegistry::add<ShapeCanvas>();
    ItemRegistry::add<Shape>();
    ItemRegistry::add<CenterShape>();
    ItemRegistry::add<RectangleShape>();
    ItemRegistry::add<PolygonShape>();
    ItemRegistry::add<ShapeDecorator>();
    ItemRegistry::add<WidthDecorator>();
    ItemRegistry::add<ColorDecorator>();
    ItemRegistry::add<ShapeGroup>();

    ItemRegistry::addUnitDef<ShapesUnitDef>();
}

////////////////////////////////////////////////////////////////////

const SchemaDef* ShapeRegistry::getSchemaDef()
{
    return ItemRegistry::getSchemaDef("Shapes");
}

////////////////////////////////////////////////////////////////////

ShapeRegistry::~ShapeRegistry()
{
    ItemRegistry::finalize();
}

////////////////////////////////////////////////////////////////////
