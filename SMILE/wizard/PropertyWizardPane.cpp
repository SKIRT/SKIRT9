/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PropertyWizardPane.hpp"
#include "PropertyHandler.hpp"
#include <QLabel>

////////////////////////////////////////////////////////////////////

PropertyWizardPane::PropertyWizardPane(std::unique_ptr<PropertyHandler> handler, QObject* target)
    : _handler(std::move(handler))
{
    connect(this, SIGNAL(propertyValidChanged(bool)), target, SLOT(setPropertyValid(bool)));
    connect(this, SIGNAL(propertyValueChanged()), target, SLOT(hierarchyWasChanged()));

    // rebuild the name sets so that conditional expressions are correctly interpreted
    _handler->rebuildNames();
}

////////////////////////////////////////////////////////////////////

PropertyWizardPane::~PropertyWizardPane() {}

////////////////////////////////////////////////////////////////////

void PropertyWizardPane::updateInterface() {}

////////////////////////////////////////////////////////////////////

void PropertyWizardPane::showEvent(QShowEvent* event)
{
    setFocus();
    focusNextChild();
    QWidget::showEvent(event);
}

////////////////////////////////////////////////////////////////////

QLabel* PropertyWizardPane::createHeader(string text)
{
    auto label = new QLabel(QString::fromStdString(text));
    label->setStatusTip(QString::fromStdString(_handler->type() + " : " + _handler->name()));
    return label;
}

////////////////////////////////////////////////////////////////////
