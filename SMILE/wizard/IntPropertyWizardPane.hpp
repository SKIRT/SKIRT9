/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef INTPROPERTYWIZARDPANE_HPP
#define INTPROPERTYWIZARDPANE_HPP

#include "PropertyWizardPane.hpp"

////////////////////////////////////////////////////////////////////

/** An IntPropertyWizardPane instance displays the user interface corresponding to an
    IntPropertyHandler. When the user changes the edit field, the corresponding value is updated in
    the target property. */
class IntPropertyWizardPane : public PropertyWizardPane
{
    Q_OBJECT

    // ============= Construction and Destruction =============

public:
    /** The default (and only) constructor creates and initializes the GUI for this pane. For a
        description of the arguments, see the PropertyWizardPane constructor. */
    explicit IntPropertyWizardPane(std::unique_ptr<PropertyHandler> handler, QObject* target);

    // ==================== Event Handling ====================

public slots:
    /** This function stores the value corresponding to the specified text string into the target
        property. */
    void updateValue(const QString& text);
};

////////////////////////////////////////////////////////////////////

#endif
