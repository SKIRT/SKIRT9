/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ITEMPROPERTYWIZARDPANE_HPP
#define ITEMPROPERTYWIZARDPANE_HPP

#include "PropertyWizardPane.hpp"
class QAbstractButton;

////////////////////////////////////////////////////////////////////

/** An ItemPropertyWizardPane instance displays the user interface corresponding to an
    ItemPropertyHandler. When the user makes a choice, a new simulation item of the selected type
    is created and stored in the target property. */
class ItemPropertyWizardPane : public PropertyWizardPane
{
    Q_OBJECT

    // ============= Construction and Destruction =============

public:
    /** The default (and only) constructor creates and initializes the GUI for this pane. For a
        description of the arguments, see the PropertyWizardPane constructor. */
    explicit ItemPropertyWizardPane(std::unique_ptr<PropertyHandler> handler, QObject* target);

    // ==================== Event Handling ====================

public slots:
    /** This function creates a simulation item of the type selected by the user and stores it into
        the target property. */
    void selectTypeFor(QAbstractButton* button);
};

////////////////////////////////////////////////////////////////////

#endif
