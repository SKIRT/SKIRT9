/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ColorDecorator.hpp"

////////////////////////////////////////////////////////////////////

void ColorDecorator::paintSelf()
{
    switch (_color)
    {
        case Color::White: setColor(1, 1, 1); break;
        case Color::Red: setColor(1, 0, 0); break;
        case Color::Green: setColor(0, 1, 0); break;
        case Color::Blue: setColor(0, 0, 1); break;
        case Color::Cyan: setColor(0, 1, 1); break;
        case Color::Magenta: setColor(1, 0, 1); break;
        case Color::Yellow: setColor(1, 1, 0); break;
        case Color::Black: setColor(0, 0, 0); break;
        case Color::Custom:
            if (_rgb.size() == 3) setColor(_rgb[0], _rgb[1], _rgb[2]);
            break;
    }
}

////////////////////////////////////////////////////////////////////
