/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DOUBLEPROPERTYWIZARDPANE_HPP
#define DOUBLEPROPERTYWIZARDPANE_HPP

#include "PropertyWizardPane.hpp"
class QLabel;
class QLineEdit;

////////////////////////////////////////////////////////////////////

/** A DoublePropertyWizardPane instance displays the user interface corresponding to a
    DoublePropertyHandler. When the user changes the edit field, the corresponding value is updated
    in the target property. */
class DoublePropertyWizardPane : public PropertyWizardPane
{
    Q_OBJECT

    // ============= Construction and Destruction =============

public:
    /** The default (and only) constructor creates and initializes the GUI for this pane. For a
        description of the arguments, see the PropertyWizardPane constructor. */
    explicit DoublePropertyWizardPane(std::unique_ptr<PropertyHandler> handler, QObject* target);

    // ==================== Event Handling ====================

protected:
    /** This function updates the interface of the wizard pane if the quantity of the handled
        property (and thus its units) have changed since construction or since the previous
        invocation of this function. */
    void updateInterface() override;

public slots:
    /** This function stores the value corresponding to the specified text string into the target
        property. */
    void updateValue(const QString& text);

    // ==================== Data Members ====================

private:
    QLabel* _header;
    QLineEdit* _field;
    string _quantity{"*"};   // the quantity currently in use, initialized to an impossible value
};

////////////////////////////////////////////////////////////////////

#endif
