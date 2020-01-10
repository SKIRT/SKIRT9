/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ShapeCanvas.hpp"
#include "Canvas.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // structure holding graphics state properties: color, width
    struct GraphicsState
    {
        double r{1.};
        double g{1.};
        double b{1.};
        double w{0.01};
    };
}

// class holding a canvas plus graphics state information
class ShapeCanvas::GraphicsStateCanvas
{
    Canvas _canvas;                 // regular canvas
    GraphicsState _cgs;             // current graphics state
    vector<GraphicsState> _gstack;  // graphics state stack
public:
    GraphicsStateCanvas(size_t numPixels) : _canvas(numPixels) {}
    void pushState() { _gstack.push_back(_cgs); }
    void popState()
    {
        _cgs = _gstack.back();
        _gstack.pop_back();
        _canvas.setColor(_cgs.r, _cgs.g, _cgs.b);
        _canvas.setWidth(_cgs.w);
    }
    void setColor(double r, double g, double b)
    {
        _canvas.setColor(r, g, b);
        _cgs.r = r;
        _cgs.g = g;
        _cgs.b = b;
    }
    void setWidth(double w)
    {
        _canvas.setWidth(w);
        _cgs.w = w;
    }
    void drawLine(double x1, double y1, double x2, double y2) { _canvas.drawLine(x1, y1, x2, y2); }
    void saveToTiff(string outPath) { _canvas.saveToTiff(outPath); }
};

////////////////////////////////////////////////////////////////////

void ShapeCanvas::paintAndSave(string outPath)
{
    try
    {
        _canvas = new GraphicsStateCanvas(500);
        _shape->paint();
        _canvas->saveToTiff(!_savePath.empty() ? _savePath : outPath);
        delete _canvas;
        _canvas = nullptr;
    }
    // ensure that canvas is deleted and pointer is cleared in all cases;
    // std::unique_ptr does not work because we need access to the canvas from subclasses
    // and we want to hide the implementation of the GraphicsStateCanvas class from the ShapeCanvas.hpp header
    catch (...)
    {
        delete _canvas;
        _canvas = nullptr;
        throw;
    }
}

////////////////////////////////////////////////////////////////////

void ShapeCanvas::pushState()
{
    if (_canvas) _canvas->pushState();
}

////////////////////////////////////////////////////////////////////

void ShapeCanvas::popState()
{
    if (_canvas) _canvas->popState();
}

////////////////////////////////////////////////////////////////////

void ShapeCanvas::setColor(double r, double g, double b)
{
    if (_canvas) _canvas->setColor(r, g, b);
}

////////////////////////////////////////////////////////////////////

void ShapeCanvas::setWidth(double w)
{
    if (_canvas) _canvas->setWidth(w);
}

////////////////////////////////////////////////////////////////////

void ShapeCanvas::drawLine(double x1, double y1, double x2, double y2)
{
    if (_canvas) _canvas->drawLine(x1, y1, x2, y2);
}

////////////////////////////////////////////////////////////////////
