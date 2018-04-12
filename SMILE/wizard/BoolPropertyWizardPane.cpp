/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BoolPropertyWizardPane.hpp"
#include "BoolPropertyHandler.hpp"
#include <QButtonGroup>
#include <QLabel>
#include <QRadioButton>
#include <QVBoxLayout>

////////////////////////////////////////////////////////////////////

BoolPropertyWizardPane::BoolPropertyWizardPane(std::unique_ptr<PropertyHandler> handler, QObject* target)
    : PropertyWizardPane(std::move(handler), target)
{
    auto hdlr = handlerCast<BoolPropertyHandler>();

    // if the property has not been configured, and there is a default value, update the property value
    if (!hdlr->isConfigured() && hdlr->hasDefaultValue())
    {
        hdlr->setValue(hdlr->defaultValue());
        hdlr->setConfigured();
    }

    // create the layout so that we can add stuff one by one
    auto layout = new QVBoxLayout;

    // add the question
    layout->addWidget(createHeader("Do you want to " + hdlr->title() + "?"));

    // make a button group to contain the radio buttons reflecting the possible choices
    auto buttonGroup = new QButtonGroup;

    // add the "No" choice
    {
        bool isDefault = hdlr->hasDefaultValue() && hdlr->defaultValue()==false;
        bool isCurrent = hdlr->isConfigured() && hdlr->value()==false;

        // add the choice button the group and to the layout
        auto choiceButton = new QRadioButton(QString("No") + (isDefault ? "  [default]" : ""));
        choiceButton->setFocusPolicy(Qt::NoFocus);
        buttonGroup->addButton(choiceButton, 0);
        layout->addWidget(choiceButton);

        // if this button corresponds to the current value, select it
        if (isCurrent)
        {
            choiceButton->setChecked(true);
            emit propertyValidChanged(true);
        }
    }

    // add the "Yes" choice
    {
        bool isDefault = hdlr->hasDefaultValue() && hdlr->defaultValue()==true;
        bool isCurrent = hdlr->isConfigured() && hdlr->value()==true;

        // add the choice button the group and to the layout
        auto choiceButton = new QRadioButton(QString("Yes") + (isDefault ? "  [default]" : ""));
        choiceButton->setFocusPolicy(Qt::NoFocus);
        buttonGroup->addButton(choiceButton, 1);
        layout->addWidget(choiceButton);

        // if this button corresponds to the current value, select it
        if (isCurrent)
        {
            choiceButton->setChecked(true);
            emit propertyValidChanged(true);
        }
    }

    // connect the button group to ourselves
    connect(buttonGroup, SIGNAL(buttonClicked(int)), this, SLOT(updateValueFor(int)));

    // finalize the layout and assign it to ourselves
    layout->addStretch();
    setLayout(layout);
}

////////////////////////////////////////////////////////////////////

void BoolPropertyWizardPane::updateValueFor(int buttonID)
{
    auto hdlr = handlerCast<BoolPropertyHandler>();

    bool newValue = buttonID ? true : false;
    if (!hdlr->isConfigured() || hdlr->value()!=newValue)
    {
        hdlr->setValue(newValue);
        emit propertyValueChanged();
    }
    hdlr->setConfigured();
    emit propertyValidChanged(true);
}

////////////////////////////////////////////////////////////////////
