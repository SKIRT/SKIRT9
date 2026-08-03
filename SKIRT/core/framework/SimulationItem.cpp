/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SimulationItem.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

void SimulationItem::setup()
{
    if (_setupStarted) return;
    _setupStarted = true;

    setupSelfBefore();
    for (Item* child : children())
    {
        SimulationItem* item = dynamic_cast<SimulationItem*>(child);
        if (item) item->setup();
    }
    setupSelfAfter();
}

////////////////////////////////////////////////////////////////////

void SimulationItem::setupSelfBefore() {}

////////////////////////////////////////////////////////////////////

void SimulationItem::setupSelfAfter() {}

////////////////////////////////////////////////////////////////////

string SimulationItem::typeAndName() const
{
    string result = type();
    string name = itemName();
    if (!name.empty()) result += " " + name;
    return result;
}

////////////////////////////////////////////////////////////////////

std::string SimulationItem::itemName() const
{
    return string();
}

////////////////////////////////////////////////////////////////////

Item* SimulationItem::find(bool setup, SimulationItem* castToRequestedType(Item*)) const
{
    // loop over all ancestors
    Item* ancestor = const_cast<SimulationItem*>(this);  // cast away const
    while (ancestor)
    {
        // test the ancestor
        SimulationItem* candidate = castToRequestedType(ancestor);
        if (candidate)
        {
            if (setup) candidate->setup();
            return candidate;
        }

        // test its children
        for (Item* child : ancestor->children())
        {
            SimulationItem* candidate = castToRequestedType(child);
            if (candidate)
            {
                if (setup) candidate->setup();
                return candidate;
            }
        }

        // next ancestor
        ancestor = ancestor->parent();
    }

    if (setup) throw FATALERROR("No simulation item of requested type found in hierarchy");
    return nullptr;
}

////////////////////////////////////////////////////////////////////

SimulationItem* SimulationItem::interface(int levels, bool setup, bool offersRequestedInterface(SimulationItem*)) const
{
    // always test the receiving object
    SimulationItem* candidate = const_cast<SimulationItem*>(this);  // cast away const
    if (offersRequestedInterface(candidate))
    {
        if (setup) candidate->setup();
        return candidate;
    }

    // test the requested number of ancestors
    if (levels < 0)
    {
        while (levels++)
        {
            candidate = dynamic_cast<SimulationItem*>(candidate->parent());
            if (!candidate) break;
            if (offersRequestedInterface(candidate))
            {
                if (setup) candidate->setup();
                return candidate;
            }
        }
    }

    // test the requested number of child levels, recursively
    else if (levels > 0)
    {
        for (Item* child : candidate->children())
        {
            auto result = dynamic_cast<SimulationItem*>(child)->interface(levels - 1, false, offersRequestedInterface);
            if (result)
            {
                if (setup) result->setup();
                return result;
            }
        }
    }

    if (setup) throw FATALERROR("No simulation item implementing requested interface found in hierarchy");
    return nullptr;
}

////////////////////////////////////////////////////////////////////

bool SimulationItem::offersInterface(const std::type_info& /*interfaceTypeInfo*/) const
{
    return true;
}

////////////////////////////////////////////////////////////////////
