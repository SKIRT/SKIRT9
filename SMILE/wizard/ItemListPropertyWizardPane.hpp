/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ITEMLISTPROPERTYWIZARDPANE_HPP
#define ITEMLISTPROPERTYWIZARDPANE_HPP

#include "PropertyWizardPane.hpp"
class QAbstractButton;
class QListWidget;
class QPushButton;

////////////////////////////////////////////////////////////////////

/** An ItemListPropertyWizardPane instance displays the user interface corresponding to an
    ItemListPropertyHandler. When the user makes a choice, a new simulation item of the selected
    type is created and added to the target property. */
class ItemListPropertyWizardPane : public PropertyWizardPane
{
    Q_OBJECT

    // ============= Construction and Destruction =============

public:
    /** The default (and only) constructor creates and initializes the GUI for this pane. For a
        description of the arguments, see the PropertyWizardPane constructor. */
    explicit ItemListPropertyWizardPane(std::unique_ptr<PropertyHandler> handler, QObject* target);

    // ==================== Event Handling ====================

public slots:
    /** This function adds a simulation item to the list of the target property. */
    void addItem();

    /** This function edits a simulation item in the list of the target property. */
    void editItem();

    /** This function removes a simulation item from the list of the target property. */
    void removeItem();

    /** This function stores the specified value for this combination of target simulation item and
        item list property. It is intended to make the zero-based index of the currently selected
        row persistent during a particular session from one display of the wizard to the next. */
    void storeSelectedRow(int row);

signals:
    /** This signal is emitted when the "Edit" button is pressed. */
    void advanceToEditSubItem(int subItemIndex);

private:
    /** This function enables or disables the push buttons depending on the contents of the list
        widget, and emits a propertyValidChanged signal with the appropriate argument. */
    void setButtonsEnabled();

    // ==================== Data members ======================

private:
    QListWidget* _listWidget;  // the widget holding the items that represent the contents of this property
    QPushButton* _addButton;
    QPushButton* _editButton;
    QPushButton* _removeButton;
};

////////////////////////////////////////////////////////////////////

#endif
