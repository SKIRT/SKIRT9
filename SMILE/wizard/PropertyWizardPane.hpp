/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PROPERTYWIZARDPANE_HPP
#define PROPERTYWIZARDPANE_HPP

#include "Basics.hpp"
#include <QWidget>
class PropertyHandler;

////////////////////////////////////////////////////////////////////

/** PropertyWizardPane is the base class for all wizard panes handling item properties. It offers
    common functionality, such as connecting the propertyValidChanged() signal to the target
    object, and retaining a reference to the relevant property handler so that it does not get
    deleted until the wizard pane is destroyed. */
class PropertyWizardPane : public QWidget
{
    Q_OBJECT

    // ============= Construction and Destruction =============

public:
    /** The default (and only) constructor retains a reference to the specified property handler so
        that it does not get deleted until the wizard pane is destroyed, and connects the
        propertyValidChanged() and propertyValueChanged() signals to the specified target object.
        */
    explicit PropertyWizardPane(std::unique_ptr<PropertyHandler> handler, QObject* target);

    /** The destructor destroys the property handler being held by this wizard pane. */
    ~PropertyWizardPane();

    // ==================== Event Handling ====================

protected:
    /** This function ensures that the first focus-enabled widget in the pane receives the focus
        when the pane is shown. */
    void showEvent(QShowEvent* event) override;

signals:
    /** This signal is emitted when the valid state of the property wizard pane (indicating whether
        the wizard is allowed to advance) may have changed. */
    void propertyValidChanged(bool valid);

    /** This signal is emitted when the value of the property being handled by this property wizard
        pane has changed. The signal should only be emitted when the value actually did change, not
        just when it might have changed. Thus the caller should check the previous property value.
        */
    void propertyValueChanged();

    // ============== Services to subclasses ==========

protected:
    /** This template function dynamically casts the handler retained by this wizard pane to the
        specified template type, and returns the result. If the handler is not of the specified
        type, the function returns null. */
    template <class T> T* handlerCast() { return dynamic_cast<T*>(_handler.get()); }

    /** This function creates a QLabel widget with the given text and with a status tip that is
        appropriate for the property being handled by this property wizard pane. */
    QWidget* createHeader(string text);

    // ================== Data Members ====================

private:
    std::unique_ptr<PropertyHandler> _handler;
};

////////////////////////////////////////////////////////////////////

#endif
