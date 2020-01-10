/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DOUBLELISTPROPERTYWIZARDPANE_HPP
#define DOUBLELISTPROPERTYWIZARDPANE_HPP

#include "PropertyWizardPane.hpp"
class QLabel;
class QLineEdit;

////////////////////////////////////////////////////////////////////

/** A DoubleListPropertyWizardPane instance displays the user interface corresponding to a
    DoublePropertyHandler. When the user changes the edit field, the corresponding value is updated
    in the target property. */
class DoubleListPropertyWizardPane : public PropertyWizardPane
{
    Q_OBJECT

    // ============= Construction and Destruction =============

public:
    /** The default (and only) constructor creates and initializes the GUI for this pane. For a
        description of the arguments, see the PropertyWizardPane constructor. */
    explicit DoubleListPropertyWizardPane(std::unique_ptr<PropertyHandler> handler, QObject* target);

    // ==================== Event Handling ====================

protected:
    /** This function updates the user interface of the pane if needed to adjust to changes to the
        values of other properties displayed inside the same MultiPropertyWizardPane instance.
        Specifically, if the quantity units and/or the minimum, maximum or default values have
        changed, the header message is updated and the field value is updated if the quantity units
        and/or the default value have changed and the field has not yet been edited by the user
        (i.e. it still contains a programmatically inserted default value). */
    void updateInterface() override;

public slots:
    /** This function stores the value corresponding to the specified text string into the target
        property. */
    void updateValue(const QString& text);

    // ==================== Data Members ====================

private:
    string _message;  // the current header message, remembered to detect changes in the environment
    QLabel* _header;
    QLineEdit* _field;
};

////////////////////////////////////////////////////////////////////

#endif
