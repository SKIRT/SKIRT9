/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SUBITEMPROPERTYWIZARDPANE_HPP
#define SUBITEMPROPERTYWIZARDPANE_HPP

#include "PropertyWizardPane.hpp"
class QAbstractButton;

////////////////////////////////////////////////////////////////////

/** A SubItemPropertyWizardPane instance displays the user interface for editing a particular
    sub-item in an item list property. When the user makes a choice, a new simulation item of the
    selected type is created and stored in the target property. */
class SubItemPropertyWizardPane : public PropertyWizardPane
{
    Q_OBJECT

    // ============= Construction and Destruction =============

public:
    /** The default (and only) constructor creates and initializes the GUI for this pane. For a
        description of the arguments, see the PropertyWizardPane constructor. */
    explicit SubItemPropertyWizardPane(std::unique_ptr<PropertyHandler> handler, QObject* target);

private:
    /** This function returns the zero-based index of the currently selected item for this
        combination of target simulation item and item list property. */
    int selectedIndex();

    // ==================== Event Handling ====================

public slots:
    /** This function creates a simulation item of the type selected by the user and stores it into
        the target property. */
    void selectTypeFor(QAbstractButton* button);
};

////////////////////////////////////////////////////////////////////

#endif
