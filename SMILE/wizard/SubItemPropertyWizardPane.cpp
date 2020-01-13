/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SubItemPropertyWizardPane.hpp"
#include "Item.hpp"
#include "ItemListPropertyHandler.hpp"
#include "SchemaDef.hpp"
#include "StringUtils.hpp"
#include <QButtonGroup>
#include <QLabel>
#include <QRadioButton>
#include <QVBoxLayout>
#include <QVariant>

////////////////////////////////////////////////////////////////////

SubItemPropertyWizardPane::SubItemPropertyWizardPane(std::unique_ptr<PropertyHandler> handler, QObject* target)
    : PropertyWizardPane(std::move(handler), target)
{
    auto hdlr = handlerCast<ItemListPropertyHandler>();

    // create the layout so that we can add stuff one by one
    auto layout = new QVBoxLayout;

    // add the question
    layout->addWidget(createHeader("Select one of the following options for item #"
                                   + std::to_string(selectedIndex() + 1) + " in " + hdlr->title() + " list:"));

    // determine the current and default item types
    string currentType = (hdlr->value()[selectedIndex()])->type();
    string defaultType = hdlr->hasDefaultValue() ? hdlr->defaultType() : "";

    // make a button group to contain the radio buttons reflecting the possible choices
    auto buttonGroup = new QButtonGroup;

    // add the choices
    for (auto choiceType : hdlr->allowedAndDisplayedDescendants())
    {
        string choiceTitle = StringUtils::toUpperFirst(hdlr->schema()->title(choiceType));
        if (hdlr->schema()->inherits(choiceType, defaultType)) choiceTitle += "  [default]";
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

int SubItemPropertyWizardPane::selectedIndex()
{
    auto hdlr = handlerCast<ItemListPropertyHandler>();
    return hdlr->retrieveSelectedRow();
}

////////////////////////////////////////////////////////////////////

void SubItemPropertyWizardPane::selectTypeFor(QAbstractButton* button)
{
    auto hdlr = handlerCast<ItemListPropertyHandler>();

    // update the value if needed
    string newType = button->property("choiceType").toString().toStdString();
    int index = selectedIndex();
    if ((hdlr->value()[index])->type() != newType)
    {
        hdlr->removeValueAt(index);
        hdlr->insertNewItemOfType(index, newType);
        emit propertyValueChanged();
    }

    // signal the change
    emit propertyValidChanged(true);
}

////////////////////////////////////////////////////////////////////
