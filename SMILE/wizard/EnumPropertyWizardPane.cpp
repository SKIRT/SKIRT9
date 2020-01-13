/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "EnumPropertyWizardPane.hpp"
#include "EnumPropertyHandler.hpp"
#include "StringUtils.hpp"
#include <QButtonGroup>
#include <QLabel>
#include <QRadioButton>
#include <QVBoxLayout>
#include <QVariant>

////////////////////////////////////////////////////////////////////

EnumPropertyWizardPane::EnumPropertyWizardPane(std::unique_ptr<PropertyHandler> handler, QObject* target)
    : PropertyWizardPane(std::move(handler), target)
{
    auto hdlr = handlerCast<EnumPropertyHandler>();

    // if the property has not been configured, and there is a default value, update the property value
    if (!hdlr->isConfigured() && hdlr->hasDefaultValue())
    {
        hdlr->setValue(hdlr->defaultValue());
        hdlr->setConfiguredToDefault();
    }

    // determine the current and default values for use in the loop below
    string currentKey = hdlr->isConfigured() ? hdlr->value() : "";
    string defaultKey = hdlr->hasDefaultValue() ? hdlr->defaultValue() : "";

    // create the layout so that we can add stuff one by one
    auto layout = new QVBoxLayout;

    // add the question
    layout->addWidget(createHeader("Select one of the following options for " + hdlr->title() + ":"));

    // make a button group to contain the radio buttons reflecting the possible choices
    auto buttonGroup = new QButtonGroup;

    // add the choices
    auto choiceKeys = hdlr->values();
    auto choiceTitles = hdlr->titlesForValues();
    for (size_t index = 0; index != choiceKeys.size(); ++index)
    {
        auto choiceKey = choiceKeys[index];
        auto choiceTitle = StringUtils::toUpperFirst(choiceTitles[index]);
        if (choiceKey == defaultKey) choiceTitle += "  [default]";
        auto choiceButton = new QRadioButton(QString::fromStdString(choiceTitle));
        choiceButton->setFocusPolicy(Qt::NoFocus);
        buttonGroup->addButton(choiceButton);
        layout->addWidget(choiceButton);

        // associate the item type corresponding to this button with the button object
        choiceButton->setProperty("choiceKey", QString::fromStdString(choiceKey));

        // if this button corresponds to the current value, select it
        if (choiceKey == currentKey)
        {
            choiceButton->setChecked(true);
            emit propertyValidChanged(true);
        }
    }

    // connect the button group to ourselves
    connect(buttonGroup, SIGNAL(buttonClicked(QAbstractButton*)), this, SLOT(updateValueFor(QAbstractButton*)));

    // finalize the layout and assign it to ourselves
    layout->addStretch();
    setLayout(layout);
}

////////////////////////////////////////////////////////////////////

void EnumPropertyWizardPane::updateValueFor(QAbstractButton* button)
{
    auto hdlr = handlerCast<EnumPropertyHandler>();

    string newValue = button->property("choiceKey").toString().toStdString();
    if (!hdlr->isConfigured() || hdlr->value() != newValue)
    {
        hdlr->setValue(newValue);
        emit propertyValueChanged();
    }
    hdlr->setConfiguredByUser();
    emit propertyValidChanged(true);
}

////////////////////////////////////////////////////////////////////
