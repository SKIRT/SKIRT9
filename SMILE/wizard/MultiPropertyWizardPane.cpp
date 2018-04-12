/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MultiPropertyWizardPane.hpp"
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
    _lastPane = nullptr;
}

////////////////////////////////////////////////////////////////////

void MultiPropertyWizardPane::addPane(PropertyWizardPane* pane)
{
    // remove the stretch layout item from the end of the layout for the previous subpane
    if (_lastPane)
    {
        auto layout = _lastPane->layout();
        auto count = layout->count();
        if (count > 0)
        {
            auto item = layout->itemAt(count-1);
            if (item->spacerItem()) layout->removeItem(item);
        }
    }

    // if it has not yet been set, initialize the subpane state to invalid
    if (!_subPaneState.contains(pane)) _subPaneState.insert(pane, false);

    // add the subpane to our layout
    _multiLayout->addWidget(pane);

    // remember this pane as the last one added
    _lastPane = pane;
}

////////////////////////////////////////////////////////////////////

void MultiPropertyWizardPane::showEvent(QShowEvent* event)
{
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
    emit multiPropertyValueChanged();
}

////////////////////////////////////////////////////////////////////
