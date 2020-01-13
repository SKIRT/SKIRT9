/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DoublePropertyWizardPane.hpp"
#include "DoublePropertyHandler.hpp"
#include <QLabel>
#include <QLineEdit>
#include <QVBoxLayout>

////////////////////////////////////////////////////////////////////

namespace
{
    // returns true if text is a valid double, and the value is within range
    bool isValidAndInRange(DoublePropertyHandler* hdlr, string text)
    {
        if (!hdlr->isValidDouble(text)) return false;
        return hdlr->isInRange(hdlr->toDouble(text));
    }

    // returns the text of the header message for this property
    string headerMessage(DoublePropertyHandler* hdlr)
    {
        string message = "Enter " + hdlr->title() + " " + hdlr->rangeDescription();
        if (hdlr->hasDefaultValue()) message += " (" + hdlr->toString(hdlr->defaultValue()) + ")";
        message += ":";
        return message;
    }

    // sets the default value for this property, if there is one
    void setDefaultValue(DoublePropertyHandler* hdlr, QLineEdit* field)
    {
        if (hdlr->hasDefaultValue())
        {
            field->setText(QString::fromStdString(hdlr->toString(hdlr->defaultValue())));
            hdlr->setValue(hdlr->defaultValue());
            hdlr->setConfiguredToDefault();
        }
        else
        {
            field->setText("");
            hdlr->setValue(0.);
            hdlr->setNotConfigured();
        }
    }
}

////////////////////////////////////////////////////////////////////

DoublePropertyWizardPane::DoublePropertyWizardPane(std::unique_ptr<PropertyHandler> handler, QObject* target)
    : PropertyWizardPane(std::move(handler), target)
{
    auto hdlr = handlerCast<DoublePropertyHandler>();

    // create the layout so that we can add stuff one by one
    auto layout = new QVBoxLayout;

    // add the header
    _message = headerMessage(hdlr);
    _header = createHeader(_message);
    layout->addWidget(_header);

    // add the edit field
    _field = new QLineEdit;
    layout->addWidget(_field);

    // connect the field to ourselves
    connect(_field, SIGNAL(textEdited(const QString&)), this, SLOT(updateValue(const QString&)));

    // finalize the layout and assign it to ourselves
    layout->addStretch();
    setLayout(layout);

    // if the property has been configured, put its value in the text field
    // otherwise, if there is a default value, put that in both the text field and the property value
    // otherwise, just empty the text field (an invalid value that will have to be updated by the user)
    if (hdlr->isConfigured())
        _field->setText(QString::fromStdString(hdlr->toString(hdlr->value())));
    else
        setDefaultValue(hdlr, _field);

    // ensure proper validity state
    emit propertyValidChanged(isValidAndInRange(hdlr, _field->text().toStdString()));
}

////////////////////////////////////////////////////////////////////

void DoublePropertyWizardPane::updateInterface()
{
    auto hdlr = handlerCast<DoublePropertyHandler>();

    // if the header message differs, this means the min/max/default values or the quantity have changed
    string newMessage = headerMessage(hdlr);
    if (newMessage != _message)
    {
        // update the header message
        _header->setText(QString::fromStdString(newMessage));
        _message = newMessage;

        // if the property has already been configured by the user,
        // verify the field value against possibly adjusted requirements (min/max value, units);
        // otherwise, (re)insert the (possibly adjusted) default value
        if (hdlr->isConfiguredByUser())
            updateValue(_field->text());
        else
            setDefaultValue(hdlr, _field);
    }
}

////////////////////////////////////////////////////////////////////

void DoublePropertyWizardPane::updateValue(const QString& text)
{
    auto hdlr = handlerCast<DoublePropertyHandler>();

    // verify that value is valid and within range before setting it
    double value = hdlr->toDouble(text.toStdString());
    bool valid = isValidAndInRange(hdlr, text.toStdString());
    if (valid && (!hdlr->isConfigured() || value != hdlr->value()))
    {
        hdlr->setValue(value);
        emit propertyValueChanged();
    }
    hdlr->setConfiguredByUser(valid);
    emit propertyValidChanged(valid);
}

////////////////////////////////////////////////////////////////////
