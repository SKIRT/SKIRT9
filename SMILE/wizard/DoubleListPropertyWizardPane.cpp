/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DoubleListPropertyWizardPane.hpp"
#include "DoubleListPropertyHandler.hpp"
#include <QLabel>
#include <QLineEdit>
#include <QVBoxLayout>

////////////////////////////////////////////////////////////////////

namespace
{
    // returns true if text is a valid double list, and all numbers are within range
    bool isValidAndInRange(DoubleListPropertyHandler* hdlr, string text)
    {
        if (!hdlr->isValidDoubleList(text)) return false;
        return hdlr->isInRange(hdlr->toDoubleList(text));
    }
}

////////////////////////////////////////////////////////////////////

DoubleListPropertyWizardPane::DoubleListPropertyWizardPane(std::unique_ptr<PropertyHandler> handler, QObject* target)
    : PropertyWizardPane(std::move(handler), target)
{
    auto hdlr = handlerCast<DoubleListPropertyHandler>();

    // create the layout so that we can add stuff one by one
    auto layout = new QVBoxLayout;

    // add the message
    string message = "Enter " + hdlr->title() + " " + hdlr->rangeDescription();
    if (hdlr->hasDefaultValue()) message += " (" + hdlr->toString(hdlr->defaultValue()) + ")";
    message += ":";
    layout->addWidget(createHeader(message));

    // add the edit field
    auto field = new QLineEdit;
    layout->addWidget(field);

    // connect the field to ourselves
    connect(field, SIGNAL(textEdited(const QString&)), this, SLOT(updateValue(const QString&)));

    // finalize the layout and assign it to ourselves
    layout->addStretch();
    setLayout(layout);

    // if the property has been configured, put its value in the text field
    // otherwise, if there is a default value, put that in both the text field and the property value
    // otherwise, just empty the text field (an invalid value that will have to be updated by the user)
    if (hdlr->isConfigured())
    {
        field->setText(QString::fromStdString(hdlr->toString(hdlr->value())));
    }
    else
    {
        if (hdlr->hasDefaultValue())
        {
            field->setText(QString::fromStdString(hdlr->toString(hdlr->defaultValue())));
            hdlr->setValue(hdlr->defaultValue());
            hdlr->setConfigured();
        }
        else
        {
            field->setText("");
        }
    }

    // ensure proper validity state
    emit propertyValidChanged(isValidAndInRange(hdlr, field->text().toStdString()));
}

////////////////////////////////////////////////////////////////////

void DoubleListPropertyWizardPane::updateValue(const QString& text)
{
    auto hdlr = handlerCast<DoubleListPropertyHandler>();

    // verify that text is valid and all numbers are within range before setting value
    auto value = hdlr->toDoubleList(text.toStdString());
    bool valid = isValidAndInRange(hdlr, text.toStdString());
    if (valid && (!hdlr->isConfigured() || value!=hdlr->value()))
    {
        hdlr->setValue(value);
        emit propertyValueChanged();
    }
    hdlr->setConfigured(valid);
    emit propertyValidChanged(valid);
}

////////////////////////////////////////////////////////////////////
