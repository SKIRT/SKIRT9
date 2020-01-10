/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ItemPropertyWizardPane.hpp"
#include "Item.hpp"
#include "ItemPropertyHandler.hpp"
#include "SchemaDef.hpp"
#include "StringUtils.hpp"
#include <QButtonGroup>
#include <QLabel>
#include <QRadioButton>
#include <QVBoxLayout>
#include <QVariant>

////////////////////////////////////////////////////////////////////

ItemPropertyWizardPane::ItemPropertyWizardPane(std::unique_ptr<PropertyHandler> handler, QObject* target)
    : PropertyWizardPane(std::move(handler), target)
{
    auto hdlr = handlerCast<ItemPropertyHandler>();

    // if there is exactly one choice, make it the forced type
    auto choiceList = hdlr->allowedAndDisplayedDescendants();
    string forcedType = (choiceList.size() == 1 && hdlr->isRequired()) ? choiceList[0] : "";

    // determine the default type, if there is one
    string defaultType = hdlr->hasDefaultValue() ? hdlr->defaultType() : "";

    // if the property has not been configured...
    if (!hdlr->isConfigured())
    {
        // if there is a forced type, update the property value to the forced type
        if (!forcedType.empty())
        {
            hdlr->setToNewItemOfType(forcedType);
            hdlr->setConfiguredToDefault();
        }
        // otherwise, if there is a default type, update the property value to the default type
        else if (!defaultType.empty())
        {
            hdlr->setToNewItemOfType(defaultType);
            hdlr->setConfiguredToDefault();
        }
        // otherwise, if the property is optional, set the property value to null
        else if (!hdlr->isRequired())
        {
            hdlr->setToNull();
            hdlr->setConfiguredToDefault();
        }
        // otherwise, leave things alone (the user will be required to select a new value anyway)
    }

    // determine the current type
    string currentType = (hdlr->isConfigured() && hdlr->value()) ? hdlr->value()->type() : "";

    // create the layout so that we can add stuff one by one
    auto layout = new QVBoxLayout;

    // add the question
    layout->addWidget(createHeader("Select one of the following options for " + hdlr->title() + ":"));

    // make a button group to contain the radio buttons reflecting the possible choices
    auto buttonGroup = new QButtonGroup;

    // if the property is optional, add the "None" choice
    if (!hdlr->isRequired())
    {
        auto choiceButton = new QRadioButton("None");
        choiceButton->setFocusPolicy(Qt::NoFocus);
        buttonGroup->addButton(choiceButton);
        layout->addWidget(choiceButton);

        // if the property currently holds no item, select this button
        if (currentType.empty())
        {
            choiceButton->setChecked(true);
            emit propertyValidChanged(true);
        }
    }

    // add the regular choices
    for (auto choiceType : choiceList)
    {
        string choiceTitle = StringUtils::toUpperFirst(hdlr->schema()->title(choiceType));
        if (choiceType == defaultType) choiceTitle += "  [default]";
        auto choiceButton = new QRadioButton(QString::fromStdString(choiceTitle));
        choiceButton->setFocusPolicy(Qt::NoFocus);
        buttonGroup->addButton(choiceButton);
        layout->addWidget(choiceButton);

        // associate the item type corresponding to this button with the button object
        choiceButton->setProperty("choiceType", QString::fromStdString(choiceType));
        choiceButton->setStatusTip(QString::fromStdString(choiceType));

        // if this button corresponds to the current type, select it
        if (choiceType == currentType)
        {
            choiceButton->setChecked(true);
            emit propertyValidChanged(true);
        }
    }

    // connect the button group to ourselves
    connect(buttonGroup, SIGNAL(buttonClicked(QAbstractButton*)), this, SLOT(selectTypeFor(QAbstractButton*)));

    // finalize the layout and assign it to ourselves
    layout->addStretch();
    setLayout(layout);
}

////////////////////////////////////////////////////////////////////

void ItemPropertyWizardPane::selectTypeFor(QAbstractButton* button)
{
    auto hdlr = handlerCast<ItemPropertyHandler>();

    // update the value
    if (button->property("choiceType").isValid())
    {
        string newType = button->property("choiceType").toString().toStdString();
        if (!hdlr->isConfigured() || !hdlr->value() || hdlr->value()->type() != newType)
        {
            hdlr->setToNewItemOfType(newType);
            emit propertyValueChanged();
        }
    }
    else
    {
        if (!hdlr->isConfigured() || hdlr->value())
        {
            hdlr->setToNull();
            emit propertyValueChanged();
        }
    }

    // make the target item remember that this property was configured by the user
    hdlr->setConfiguredByUser();

    // signal the change
    emit propertyValidChanged(true);
}

////////////////////////////////////////////////////////////////////
