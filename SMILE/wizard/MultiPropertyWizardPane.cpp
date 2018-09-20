/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MultiPropertyWizardPane.hpp"
#include "PropertyHandler.hpp"
#include "PropertyWizardPane.hpp"
#include <QHash>
#include <QVBoxLayout>
#include <QSpacerItem>

////////////////////////////////////////////////////////////////////

MultiPropertyWizardPane::MultiPropertyWizardPane(QObject* target)
{
    // connect our signals to the target
    connect(this, SIGNAL(multiPropertyValidChanged(bool)), target, SLOT(setPropertyValid(bool)));
    connect(this, SIGNAL(multiPropertyValueChanged()), target, SLOT(hierarchyWasChanged()));

    // create a layout and assign it to ourselves
    _multiLayout = new QVBoxLayout;
    _multiLayout->setContentsMargins(0,0,0,0);
    _multiLayout->setSpacing(0);
    setLayout(_multiLayout);
}

////////////////////////////////////////////////////////////////////

void MultiPropertyWizardPane::addPane(PropertyWizardPane* pane)
{
    // remove the stretch layout item from the end of the layout for the previous subpane
    if (!_subPanes.empty())
    {
        auto layout = _subPanes.back()->layout();
        auto count = layout->count();
        if (count > 0)
        {
            auto item = layout->itemAt(count-1);
            if (item->spacerItem()) layout->removeItem(item);
        }
    }

    // add the subpane to our layout
    _multiLayout->addWidget(pane);

    // add the subpane to our list
    _subPanes.append(pane);

    // if it has not yet been set, initialize the subpane state to invalid
    if (!_subPaneState.contains(pane)) _subPaneState.insert(pane, false);
}

////////////////////////////////////////////////////////////////////

void MultiPropertyWizardPane::updateVisibility()
{
    if (!_subPanes.empty())
    {
        // initialize name manager up to just before the property in the first subpane
        _subPanes.front()->handler()->rebuildNames();

        // loop over all subpanes
        for (auto pane : _subPanes)
        {
            auto handler = pane->handler();

            // a property is silent if it is irrelevant,
            // or if it should not be displayed (unless it is required and has no default value)
            bool silent = handler->isSilent();

            // set the visibility accordingly
            if (!silent) pane->updateInterface();
            pane->setVisible(!silent);

            // insert the names for the property
            handler->insertNames();
        }
    }
}

////////////////////////////////////////////////////////////////////

void MultiPropertyWizardPane::showEvent(QShowEvent* event)
{
    updateVisibility();
    setFocus();
    focusNextChild();
    QWidget::showEvent(event);
}

////////////////////////////////////////////////////////////////////

void MultiPropertyWizardPane::setPropertyValid(bool valid)
{
    // update subpane state
    _subPaneState.insert(sender(), valid);

    // calculate and emit the aggregate state
    auto states = _subPaneState.values();
    bool multivalid = std::all_of(states.cbegin(), states.cend(), [](bool s){return s;});
    emit multiPropertyValidChanged(multivalid);
}

////////////////////////////////////////////////////////////////////

void MultiPropertyWizardPane::hierarchyWasChanged()
{
    updateVisibility();
    emit multiPropertyValueChanged();
}

////////////////////////////////////////////////////////////////////
