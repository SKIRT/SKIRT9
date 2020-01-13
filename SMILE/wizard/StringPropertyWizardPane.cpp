/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StringPropertyWizardPane.hpp"
#include "StringPropertyHandler.hpp"
#include <QLabel>
#include <QLineEdit>
#include <QVBoxLayout>

////////////////////////////////////////////////////////////////////

StringPropertyWizardPane::StringPropertyWizardPane(std::unique_ptr<PropertyHandler> handler, QObject* target)
    : PropertyWizardPane(std::move(handler), target)
{
    auto hdlr = handlerCast<StringPropertyHandler>();

    // create the layout so that we can add stuff one by one
    auto layout = new QVBoxLayout;

    // add the message
    string message = "Enter " + hdlr->title();
    if (hdlr->hasDefaultValue()) message += " (" + hdlr->defaultValue() + ")";
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
    // otherwise, if the property is optional, put an empty string in both the text field and the property value
    // otherwise, just empty the text field (an invalid value that will have to be updated by the user)
    if (hdlr->isConfigured())
    {
        field->setText(QString::fromStdString(hdlr->value()));
    }
    else
    {
        if (hdlr->hasDefaultValue() || !hdlr->isRequired())
        {
            field->setText(QString::fromStdString(hdlr->defaultValue()));
            hdlr->setValue(hdlr->defaultValue());
            hdlr->setConfiguredToDefault();
        }
        else
        {
            field->setText("");
        }
    }

    // ensure proper validity state
    emit propertyValidChanged(!hdlr->isRequired() || !field->text().isEmpty());
}

////////////////////////////////////////////////////////////////////

void StringPropertyWizardPane::updateValue(const QString& text)
{
    auto hdlr = handlerCast<StringPropertyHandler>();

    // verify that value is non-empty before setting it
    string value = text.simplified().toStdString();
    bool valid = !hdlr->isRequired() || !value.empty();
    if (valid && (!hdlr->isConfigured() || value != hdlr->value()))
    {
        hdlr->setValue(value);
        emit propertyValueChanged();
    }
    hdlr->setConfiguredByUser(valid);
    emit propertyValidChanged(valid);
}

////////////////////////////////////////////////////////////////////
