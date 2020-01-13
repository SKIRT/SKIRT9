/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ItemListPropertyWizardPane.hpp"
#include "Item.hpp"
#include "ItemListPropertyHandler.hpp"
#include "ItemPropertyHandler.hpp"
#include "ItemUtils.hpp"
#include "SchemaDef.hpp"
#include "StringPropertyHandler.hpp"
#include <QHBoxLayout>
#include <QLabel>
#include <QListWidget>
#include <QPushButton>
#include <QVBoxLayout>
#include <QVariant>

////////////////////////////////////////////////////////////////////

namespace
{
    // return a description for the specified simulation item containing the item type
    // and the name and value of the first item or string property in the item
    string descriptionForItem(const SchemaDef* schema, Item* item)
    {
        for (const string& property : schema->properties(item->type()))
        {
            auto handler = schema->createPropertyHandler(item, property, nullptr);
            auto itemhandler = dynamic_cast<ItemPropertyHandler*>(handler.get());
            if (itemhandler && itemhandler->isConfigured() && itemhandler->value())
            {
                return item->type() + " (" + itemhandler->name() + ": " + itemhandler->value()->type() + ")";
            }
            auto stringhandler = dynamic_cast<StringPropertyHandler*>(handler.get());
            if (stringhandler && stringhandler->isConfigured() && !stringhandler->value().empty())
            {
                return item->type() + " (" + stringhandler->name() + ": " + stringhandler->value() + ")";
            }
        }
        return item->type();
    }
}

////////////////////////////////////////////////////////////////////

ItemListPropertyWizardPane::ItemListPropertyWizardPane(std::unique_ptr<PropertyHandler> handler, QObject* target)
    : PropertyWizardPane(std::move(handler), target)
{
    auto hdlr = handlerCast<ItemListPropertyHandler>();

    // connect our edit-sub-item signal to the wizard
    connect(this, SIGNAL(advanceToEditSubItem(int)), target, SLOT(advanceToEditSubItem(int)));

    // ---- construct the overall layout ----

    // create the layout so that we can add stuff one by one
    auto layout = new QVBoxLayout;

    // add the question
    layout->addWidget(createHeader("Compose " + hdlr->title() + " list:"));

    // add the widget holding the items that represent the contents of this property
    _listWidget = new QListWidget;
    _listWidget->setSelectionMode(QAbstractItemView::SingleSelection);
    layout->addWidget(_listWidget);

    // connect the widget to our respective slot
    connect(_listWidget, SIGNAL(currentRowChanged(int)), this, SLOT(storeSelectedRow(int)));

    // add the push buttons for adding, editing and removing items
    _removeButton = new QPushButton("Remove");
    _editButton = new QPushButton("Edit");
    _addButton = new QPushButton("Add");
    auto buttonLayout = new QHBoxLayout;
    layout->addLayout(buttonLayout);
    buttonLayout->addWidget(_removeButton);
    buttonLayout->addWidget(_editButton);
    buttonLayout->addWidget(_addButton);

    // connect the buttons to our respective slots
    connect(_addButton, SIGNAL(clicked()), this, SLOT(addItem()));
    connect(_editButton, SIGNAL(clicked()), this, SLOT(editItem()));
    connect(_removeButton, SIGNAL(clicked()), this, SLOT(removeItem()));

    // finalize the layout and assign it to ourselves
    setLayout(layout);

    // ---- populate the list widget with its initial contents ----

    // if the property has not been configured, set the property value to an empty list
    if (!hdlr->isConfigured())
    {
        hdlr->setToEmpty();
        hdlr->setConfiguredToDefault();
    }

    // iterate over all items in the item list
    int index = 0;
    for (auto item : hdlr->value())
    {
        index++;
        string label = std::to_string(index) + ": " + descriptionForItem(hdlr->schema(), item);
        if (!ItemUtils::isItemComplete(item)) label += u8"  \u2190 editing incomplete !";
        _listWidget->addItem(new QListWidgetItem(QString::fromStdString(label)));
    }

    // restore the selection
    int row = hdlr->retrieveSelectedRow();
    if (row < _listWidget->count()) _listWidget->setCurrentRow(row);

    // enable buttons appropriately
    setButtonsEnabled();
}

////////////////////////////////////////////////////////////////////

void ItemListPropertyWizardPane::addItem()
{
    auto hdlr = handlerCast<ItemListPropertyHandler>();

    // select a new item type that inherits from the default, or if none is found, use the first type in the list
    auto choiceList = hdlr->allowedAndDisplayedDescendants();
    if (choiceList.empty()) return;  // can we do something more informative to the user in case of failure?
    auto newType = choiceList[0];
    if (hdlr->hasDefaultValue())
    {
        auto defaultType = hdlr->defaultType();
        for (auto choiceType : choiceList)
        {
            if (hdlr->schema()->inherits(choiceType, defaultType))
            {
                newType = choiceType;
                break;
            }
        }
    }

    // add a new item of the selected type to the property's list
    bool success = hdlr->addNewItemOfType(newType);
    if (!success) return;  // can we do something more informative to the user in case of failure?
    emit propertyValueChanged();

    // add a corresponding line to the list widget
    int count = _listWidget->count();
    _listWidget->addItem(new QListWidgetItem(QString::fromStdString(std::to_string(count + 1) + ": " + newType)));
    _listWidget->setCurrentRow(count);

    // start the item edit wizard for the current row
    editItem();
}

////////////////////////////////////////////////////////////////////

void ItemListPropertyWizardPane::editItem()
{
    auto hdlr = handlerCast<ItemListPropertyHandler>();
    emit advanceToEditSubItem(hdlr->retrieveSelectedRow());
}

////////////////////////////////////////////////////////////////////

void ItemListPropertyWizardPane::removeItem()
{
    if (_listWidget->count() > 0)
    {
        // remove from widget
        int index = _listWidget->currentRow();
        delete _listWidget->takeItem(index);

        // remove from item
        auto hdlr = handlerCast<ItemListPropertyHandler>();
        hdlr->removeValueAt(index);

        // explicitly store the newly selected row index because the widget signal sometimes messes up
        storeSelectedRow(_listWidget->currentRow());

        // tell the world that the property value changed
        emit propertyValueChanged();
    }
    setButtonsEnabled();
}

////////////////////////////////////////////////////////////////////

void ItemListPropertyWizardPane::storeSelectedRow(int row)
{
    auto hdlr = handlerCast<ItemListPropertyHandler>();
    hdlr->storeSelectedRow(row);
}

////////////////////////////////////////////////////////////////////

void ItemListPropertyWizardPane::setButtonsEnabled()
{
    // we assume that an item is always selected unless the list is empty
    bool hasItems = _listWidget->count() > 0;

    // check whether all items are completed
    auto hdlr = handlerCast<ItemListPropertyHandler>();
    bool complete = true;
    for (auto item : hdlr->value())
        if (!ItemUtils::isItemComplete(item)) complete = false;

    // enable/disable buttons
    _removeButton->setEnabled(hasItems);
    _editButton->setEnabled(hasItems);
    _addButton->setEnabled(complete);  // block new editing because the name manager dislikes incomplete items

    // emit validate/invalidate signal
    emit propertyValidChanged(complete && (hasItems || !hdlr->isRequired()));
}

////////////////////////////////////////////////////////////////////
