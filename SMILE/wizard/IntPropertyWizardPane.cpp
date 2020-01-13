/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "IntPropertyWizardPane.hpp"
#include "IntPropertyHandler.hpp"
#include "StringUtils.hpp"
#include <QLabel>
#include <QLineEdit>
#include <QVBoxLayout>

////////////////////////////////////////////////////////////////////

namespace
{
    // returns true if text is a valid integer, and the value is within range
    bool isValidAndInRange(IntPropertyHandler* hdlr, string text)
    {
        if (!StringUtils::isValidInt(text)) return false;
        int value = StringUtils::toInt(text);
        if (value < hdlr->minValue() || value > hdlr->maxValue()) return false;
        return true;
    }
}

////////////////////////////////////////////////////////////////////

IntPropertyWizardPane::IntPropertyWizardPane(std::unique_ptr<PropertyHandler> handler, QObject* target)
    : PropertyWizardPane(std::move(handler), target)
{
    auto hdlr = handlerCast<IntPropertyHandler>();

    // create the layout so that we can add stuff one by one
    auto layout = new QVBoxLayout;

    // add the message
    string message = "Enter " + hdlr->title();
    message += " [" + StringUtils::toString(hdlr->minValue()) + "," + StringUtils::toString(hdlr->maxValue()) + "]";
    if (hdlr->hasDefaultValue()) message += " (" + StringUtils::toString(hdlr->defaultValue()) + ")";
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
        field->setText(QString::fromStdString(StringUtils::toString(hdlr->value())));
    }
    else
    {
        if (hdlr->hasDefaultValue())
        {
            field->setText(QString::fromStdString(StringUtils::toString(hdlr->defaultValue())));
            hdlr->setValue(hdlr->defaultValue());
            hdlr->setConfiguredToDefault();
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

void IntPropertyWizardPane::updateValue(const QString& text)
{
    auto hdlr = handlerCast<IntPropertyHandler>();

    // verify that value is valid and within range before setting it
    int value = StringUtils::toInt(text.toStdString());
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
