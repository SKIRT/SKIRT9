/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STRINGPROPERTYWIZARDPANE_HPP
#define STRINGPROPERTYWIZARDPANE_HPP

#include "PropertyWizardPane.hpp"

////////////////////////////////////////////////////////////////////

/** A StringPropertyWizardPane instance displays the user interface corresponding to a
    StringPropertyHandler. When the user changes the edit field, the corresponding value is updated
    in the target property. */
class StringPropertyWizardPane : public PropertyWizardPane
{
    Q_OBJECT

    // ============= Construction and Destruction =============

public:
    /** The default (and only) constructor creates and initializes the GUI for this pane. For a
        description of the arguments, see the PropertyWizardPane constructor. */
    explicit StringPropertyWizardPane(std::unique_ptr<PropertyHandler> handler, QObject* target);

    // ==================== Event Handling ====================

public slots:
    /** This function stores the value corresponding to the specified text string into the target
        property. */
    void updateValue(const QString& text);
};

////////////////////////////////////////////////////////////////////

#endif
