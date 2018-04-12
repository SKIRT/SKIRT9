/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "WizardEngine.hpp"
#include "BasicChoiceWizardPane.hpp"
#include "BoolPropertyHandler.hpp"
#include "BoolPropertyWizardPane.hpp"
#include "CreateRootWizardPane.hpp"
#include "DoubleListPropertyHandler.hpp"
#include "DoubleListPropertyWizardPane.hpp"
#include "DoublePropertyHandler.hpp"
#include "DoublePropertyWizardPane.hpp"
#include "EnumPropertyHandler.hpp"
#include "EnumPropertyWizardPane.hpp"
#include "IntPropertyHandler.hpp"
#include "IntPropertyWizardPane.hpp"
#include "Item.hpp"
#include "ItemListPropertyHandler.hpp"
#include "ItemListPropertyWizardPane.hpp"
#include "ItemPropertyHandler.hpp"
#include "ItemPropertyWizardPane.hpp"
#include "ItemUtils.hpp"
#include "MultiPropertyWizardPane.hpp"
#include "OpenWizardPane.hpp"
#include "PropertyHandlerVisitor.hpp"
#include "SaveWizardPane.hpp"
#include "SchemaDef.hpp"
#include "StringPropertyHandler.hpp"
#include "StringPropertyWizardPane.hpp"
#include "StringUtils.hpp"
#include "SubItemPropertyWizardPane.hpp"

////////////////////////////////////////////////////////////////////

WizardEngine::WizardEngine(QObject* parent)
    : QObject(parent)
{
}

////////////////////////////////////////////////////////////////////

WizardEngine::~WizardEngine()
{
}

////////////////////////////////////////////////////////////////////

bool WizardEngine::canAdvance()
{
    switch (_state)
    {
    case State::BasicChoice:
        return _schema != nullptr;
    case State::CreateRoot:
        return (_schema && _root) ? _schema->inherits(_root->type(), _schema->schemaType()) : false;
    case State::OpenHierarchy:
        return !_filepath.isEmpty();
    case State::ConstructHierarchy:
        return _propertyValid;
    case State::SaveHierarchy:
        return false;
    }
    return false;   // to satisfy some compilers
}

////////////////////////////////////////////////////////////////////

bool WizardEngine::canRetreat()
{
    return _state != State::BasicChoice;
}

////////////////////////////////////////////////////////////////////

bool WizardEngine::isDirty()
{
    return _dirty;
}

////////////////////////////////////////////////////////////////////

QString WizardEngine::filepath()
{
    return _filepath;
}

////////////////////////////////////////////////////////////////////

const SchemaDef* WizardEngine::schema()
{
    return _schema.get();
}

////////////////////////////////////////////////////////////////////

int WizardEngine::propertyIndexForChild(Item* parent, Item* child)
{
    int index = 0;
    for (auto property : _schema->properties(parent->type()))
    {
        auto handler = _schema->createPropertyHandler(parent, property);

        // check the value of item properties
        auto itemhandler = dynamic_cast<ItemPropertyHandler*>(handler.get());
        if (itemhandler && itemhandler->value() == child) return index;

        // check the values of item list properties
        auto itemlisthandler = dynamic_cast<ItemListPropertyHandler*>(handler.get());
        if (itemlisthandler)
        {
            for (auto item : itemlisthandler->value()) if (item == child) return index;
        }

        index++;
    }
    return 0;  // this should never happen
}

////////////////////////////////////////////////////////////////////

std::unique_ptr<PropertyHandler> WizardEngine::createCurrentPropertyHandler(int offset)
{
    return _schema->createPropertyHandler(_current, _schema->properties(_current->type())[_propertyIndex+offset]);
}

////////////////////////////////////////////////////////////////////

namespace
{
    // Returns true if the specified property is "silent" for purposes of the wizard; this includes
    // silent properties, irrelevant properties, item properties that offer only one choice,
    // and the subitem pane of item list properties that offer only one choice (if the second argument is provided)
    bool silentForWizard(PropertyHandler* handler, int subItemIndex = -1)
    {
        if (handler->isSilent() || !handler->isRelevant()) return true;

        auto itemhdlr = dynamic_cast<ItemPropertyHandler*>(handler);
        if (itemhdlr && !itemhdlr->isOptional() && itemhdlr->allowedDescendants().size() == 1) return true;

        if (subItemIndex < 0) return false;
        auto itemlisthdlr = dynamic_cast<ItemListPropertyHandler*>(handler);
        if (itemlisthdlr && itemlisthdlr->allowedDescendants().size() == 1) return true;

        return false;
    }

    // The functions in this class are part of the visitor pattern initiated by the setupProperties() function.
    // They set the value of a non-compound property to its default value, the value of a compound property
    // to "not present", i.e. the null pointer or the empty list, and the value of an item property that
    // offers only one choice to the "forced" value
    class SilentPropertySetter : public PropertyHandlerVisitor
    {
    public:
        SilentPropertySetter() { }

        void visitPropertyHandler(StringPropertyHandler* handler) override
        {
            handler->setValue(handler->defaultValue());
        }

        void visitPropertyHandler(BoolPropertyHandler* handler) override
        {
            handler->setValue(handler->defaultValue());
        }

        void visitPropertyHandler(IntPropertyHandler* handler) override
        {
            handler->setValue(handler->defaultValue());
        }

        void visitPropertyHandler(EnumPropertyHandler* handler) override
        {
            handler->setValue(handler->defaultValue());
        }

        void visitPropertyHandler(DoublePropertyHandler* handler) override
        {
            handler->setValue(handler->defaultValue());
        }

        void visitPropertyHandler(DoubleListPropertyHandler* handler) override
        {
            handler->setValue(handler->defaultValue());
        }

        void visitPropertyHandler(ItemPropertyHandler* handler) override
        {
            if (!handler->isOptional() && handler->allowedDescendants().size() == 1)
                handler->setToNewItemOfType(handler->allowedDescendants()[0]);
            else
                handler->setToNull();
        }

        void visitPropertyHandler(ItemListPropertyHandler* handler) override
        {
            handler->setToEmpty();
        }
    };
}

////////////////////////////////////////////////////////////////////

void WizardEngine::advance()
{
    switch (_state)
    {
    case State::BasicChoice:
        {
            _state = _openExisting ? State::OpenHierarchy : State::CreateRoot;
        }
        break;
    case State::OpenHierarchy:
        {
            _state = State::CreateRoot;
        }
        break;
    case State::CreateRoot:
        {
            _state = State::ConstructHierarchy;
            _current = _root.get();
            _propertyIndex = 0; // assumes that the root has at least one property
            _propertyValid = false;
            // indicate not editing sub-item (meaningless and harmless if current item is not an item list)
            _subItemIndex = -1;
        }
        break;
    case State::ConstructHierarchy:
        {
            auto handler = createCurrentPropertyHandler();
            auto itemhdlr = dynamic_cast<ItemPropertyHandler*>(handler.get());
            auto itemlisthdlr = dynamic_cast<ItemListPropertyHandler*>(handler.get());

            // if the property being handled is an item, and the item has properties, then descend the hierarchy
            if (itemhdlr && itemhdlr->value() && _schema->properties(itemhdlr->value()->type()).size()>0)
            {
                _current = itemhdlr->value();
                _propertyIndex = 0;
            }

            // otherwise, if the property being handled is an item list, and we're editing one of its subitems,
            // and the subitem has properties, then descend the hierarchy into that subitem
            else if (itemlisthdlr && _subItemIndex>=0 &&
                     _schema->properties(itemlisthdlr->value()[_subItemIndex]->type()).size()>0)
            {
                _current = itemlisthdlr->value()[_subItemIndex];
                _propertyIndex = 0;
            }

            // otherwise, if this was the last property at this level, move up the hierarchy to a level where
            // there are properties to advance to; if we encounter the root item, then move to the SaveHierarchy state
            else
            {
                while (static_cast<size_t>(_propertyIndex+1) == _schema->properties(_current->type()).size())
                {
                    // indicate that the item we're backing out of is "complete"
                    ItemUtils::setItemComplete(_current);

                    // special case for root
                    if (_current == _root.get())
                    {
                        _state = State::SaveHierarchy;
                        break;
                    }

                    // move up the hierarchy
                    Item* parent = _current->parent();
                    _propertyIndex = propertyIndexForChild(parent, _current);
                    _current = parent;

                    // if we're advancing out of a subitem, go to the item list property rather than the next property
                    if (dynamic_cast<ItemListPropertyHandler*>(createCurrentPropertyHandler().get())) _propertyIndex--;
                }

                // advance to the next property (meaningless and harmless if state has changed to SaveHierarchy)
                _propertyIndex++;
            }

            // indicate property invalid (meaningless and harmless if state has changed to SaveHierarchy)
            _propertyValid = false;

            // indicate not editing sub-item (meaningless and harmless if current item is not an item list)
            _subItemIndex = -1;
        }
        break;
    case State::SaveHierarchy:
        break;
    }

    if (_state == State::CreateRoot)
    {
        // skip create root pane if it offers only one choice
        auto choices = _schema->descendants(_schema->schemaType());
        if (choices.size() == 1)
        {
            if (!_root || !_schema->inherits(_root->type(), _schema->schemaType()))
            {
                setRootType(choices[0]);
            }
            return advance();
        }
    }
    else if (_state == State::ConstructHierarchy)
    {
        // skip silent properties, irrelevant properties, and item properties that offer only one choice
        auto handler = createCurrentPropertyHandler();
        if (silentForWizard(handler.get()))
        {
            if (!handler->isConfigured())
            {
                SilentPropertySetter silentSetter;
                handler->acceptVisitor(&silentSetter);
                handler->setConfigured();
            }
            return advance();
        }
    }

    emitStateChanged();
}

////////////////////////////////////////////////////////////////////

void WizardEngine::descendToDeepest()
{
    while (true)
    {
        // if the current property is an item, and the item has properties, then descend the hierarchy
        auto handler = createCurrentPropertyHandler();
        auto itemhandler = dynamic_cast<ItemPropertyHandler*>(handler.get());
        if (itemhandler && itemhandler->value() && _schema->properties(itemhandler->value()->type()).size()>0)
        {
            _current = itemhandler->value();
            _propertyIndex = static_cast<int>(_schema->properties(_current->type()).size()) - 1;
        }
        else break;
    }
}

////////////////////////////////////////////////////////////////////

namespace
{
    // returns the index of a given item in the list
    template <class T> int indexOf(const vector<T>& list, const T& item)
    {
        auto pos = find(list.cbegin(), list.cend(), item);
        return pos != list.cend() ? pos - list.cbegin() : -1;
    }

    // returns true if the specified property is of type Bool, Int, or Double,
    // it is not silent and it cannot become irrelevant
    bool eligableForMultiPane(PropertyHandler* handler)
    {
        if (handler->isSilent()) return false;
        if (handler->hasRelevantIf()) return false;
        if (dynamic_cast<BoolPropertyHandler*>(handler)) return true;
        if (dynamic_cast<IntPropertyHandler*>(handler)) return true;
        if (dynamic_cast<DoublePropertyHandler*>(handler)) return true;
        return false;
    }
}

////////////////////////////////////////////////////////////////////

void WizardEngine::retreat()
{
    switch (_state)
    {
    case State::BasicChoice:
        break;
    case State::OpenHierarchy:
        {
            _state = State::BasicChoice;
        }
        break;
    case State::CreateRoot:
        {
            _state = _openExisting ? State::OpenHierarchy : State::BasicChoice;
        }
        break;
    case State::ConstructHierarchy:
        {
            // if we're retreating from a multi-pane, scuttle back to the first property on the pane
            if (eligableForMultiPane(createCurrentPropertyHandler().get()))
            {
                while (_propertyIndex>0)
                {
                    if (!eligableForMultiPane(createCurrentPropertyHandler(-1).get())) break;
                    _propertyIndex--;
                }
            }

            // if this is an item list property, and we're editing a subitem,
            // go to the item list property rather than the previous property
            if (dynamic_cast<ItemListPropertyHandler*>(createCurrentPropertyHandler().get())
                     && _subItemIndex>=0) _subItemIndex = -1;

            // otherwise, if this was the first property at this level, move up the hierarchy to the previous level
            // unless this is already the root item, in which case we move to the CreateRoot state
            else if (_propertyIndex == 0)
            {
                if (_current == _root.get())
                {
                    _state = State::CreateRoot;
                }
                else
                {
                    auto child = _current;
                    _current = child->parent();
                    _propertyIndex = propertyIndexForChild(_current, child);

                    // if we're retreating out of a subitem, go to the subitem choice pane first
                    auto handler = createCurrentPropertyHandler();
                    auto itemlisthdlr = dynamic_cast<ItemListPropertyHandler*>(handler.get());
                    if (itemlisthdlr) _subItemIndex = indexOf(itemlisthdlr->value(), child);
                    else _subItemIndex = -1;
                }
            }

            // otherwise, retreat to the previous property, and descend its hierarchy as deep as possible
            else
            {
                _propertyIndex--;
                descendToDeepest();

                // indicate not editing sub-item (meaningless and harmless if current item is not an item list)
                _subItemIndex = -1;
            }

            // indicate property invalid (meaningless and harmless if state has changed to CreateRoot)
            _propertyValid = false;
        }
        break;
    case State::SaveHierarchy:
        {
            // go back to hierarchy construction
            _state = State::ConstructHierarchy;
            _propertyValid = false;

            // descend the existing hierarchy as deep as possible
            // assumes that the root has at least one property
            _current = _root.get();
            _propertyIndex = static_cast<int>(_schema->properties(_root->type()).size()) - 1;
            descendToDeepest();

            // indicate not editing sub-item (meaningless and harmless if current item is not an item list)
            _subItemIndex = -1;
        }
        break;
    }

    if (_state == State::ConstructHierarchy)
    {
        // skip silent, irrelevant and "forced value" properties
        auto handler = createCurrentPropertyHandler();
        if (silentForWizard(handler.get(), _subItemIndex)) return retreat();

        // if we landed on a multi-pane, scuttle back to the first property on the pane
        if (eligableForMultiPane(createCurrentPropertyHandler().get()))
        {
            while (_propertyIndex>0)
            {
                if (!eligableForMultiPane(createCurrentPropertyHandler(-1).get())) break;
                _propertyIndex--;
            }
        }
    }
    else if (_state == State::CreateRoot)
    {
        // skip create root pane if it offers only one choice
        if (_schema->descendants(_schema->schemaType()).size() == 1) return retreat();
    }

    emitStateChanged();
}

////////////////////////////////////////////////////////////////////

void WizardEngine::advanceToEditSubItem(int subItemIndex)
{
    // indicate that we're editing the current sub-item
    _subItemIndex = subItemIndex;

    // skip this wizard pane if there is only one choice for the subitem class
    auto handler = createCurrentPropertyHandler();
    if (silentForWizard(handler.get(), _subItemIndex)) return advance();

    emitStateChanged();
}

////////////////////////////////////////////////////////////////////

void WizardEngine::emitStateChanged()
{
    emit stateChanged();
    emit canAdvanceChangedTo(canAdvance());
    emit canRetreatChangedTo(canRetreat());
}

////////////////////////////////////////////////////////////////////

void WizardEngine::setBasicChoice(bool openExisting, string libraryPath, string schemaName)
{
    if (_openExisting != openExisting || _schemaName != schemaName)
    {
        // update the choice
        _openExisting = openExisting;
        _schemaName = schemaName;
        _schema = std::make_unique<SchemaDef>(StringUtils::joinPaths(libraryPath, schemaName));

        // clear the current hierarchy and the related state
        _root.reset();
        _filepath.clear();
        _dirty = false;
        emit titleChanged();
        emit dirtyChanged();
        emit canAdvanceChangedTo(canAdvance());
    }
}

////////////////////////////////////////////////////////////////////

void WizardEngine::setRootType(string newRootType)
{
    if (_root && _root->type()==newRootType) return;
    _root = _schema->createItem(newRootType);
    emit canAdvanceChangedTo(canAdvance());
    _dirty = true;
    emit dirtyChanged();
}

////////////////////////////////////////////////////////////////////

void WizardEngine::hierarchyWasLoaded(Item* root, QString filepath)
{
    _root.reset(root);
    hierarchyWasSaved(filepath);
}

////////////////////////////////////////////////////////////////////

void WizardEngine::setPropertyValid(bool valid)
{
    _propertyValid = valid;
    emit canAdvanceChangedTo(canAdvance());
}

////////////////////////////////////////////////////////////////////

void WizardEngine::hierarchyWasChanged()
{
    _dirty = true;
    emit dirtyChanged();
    ItemUtils::setItemIncomplete(_current);
}

////////////////////////////////////////////////////////////////////

void WizardEngine::hierarchyWasSaved(QString filepath)
{
    _filepath = filepath;
    _dirty = false;
    emit titleChanged();
    emit dirtyChanged();
    emit canAdvanceChangedTo(canAdvance());
}

////////////////////////////////////////////////////////////////////

QWidget* WizardEngine::createPane()
{
    switch (_state)
    {
    case State::BasicChoice:
        {
            return new BasicChoiceWizardPane(_openExisting, _schemaName, _dirty, this);
        }
    case State::CreateRoot:
        {
            string currentType = _root ? _root->type() : "";
            return new CreateRootWizardPane(_schema.get(), currentType, this);
        }
    case State::OpenHierarchy:
        {
            return new OpenWizardPane(_schema.get(), _filepath, _dirty, this);
        }
    case State::ConstructHierarchy:
        {
            auto handler = createCurrentPropertyHandler();

            // handle properties that are not combined on a multi-pane
            if (dynamic_cast<StringPropertyHandler*>(handler.get()))
                return new StringPropertyWizardPane(std::move(handler), this);
            if (dynamic_cast<EnumPropertyHandler*>(handler.get()))
                return new EnumPropertyWizardPane(std::move(handler), this);
            if (dynamic_cast<DoubleListPropertyHandler*>(handler.get()))
                return new DoubleListPropertyWizardPane(std::move(handler), this);
            if (dynamic_cast<ItemPropertyHandler*>(handler.get()))
            {
                return new ItemPropertyWizardPane(std::move(handler), this);
            }
            if (dynamic_cast<ItemListPropertyHandler*>(handler.get()))
            {
                if (_subItemIndex<0) return new ItemListPropertyWizardPane(std::move(handler), this);
                else                 return new SubItemPropertyWizardPane(std::move(handler), this);
            }

            // if we reach here, we have a Bool, Int or Double property;
            // gather subsequent properties of these types at the same level in a single multi-pane
            int numProperties = static_cast<int>(_schema->properties(_current->type()).size());
            auto multipane = new MultiPropertyWizardPane(this);
            while (true)
            {
                // create and add the appropriate pane for the current property
                PropertyWizardPane* pane = nullptr;
                if (dynamic_cast<BoolPropertyHandler*>(handler.get()))
                    pane = new BoolPropertyWizardPane(std::move(handler), multipane);
                else if (dynamic_cast<IntPropertyHandler*>(handler.get()))
                    pane = new IntPropertyWizardPane(std::move(handler), multipane);
                else if (dynamic_cast<DoublePropertyHandler*>(handler.get()))
                    pane = new DoublePropertyWizardPane(std::move(handler), multipane);
                if (pane) multipane->addPane(pane);

                // terminate if this was the last property at this level
                if (_propertyIndex+1 >= numProperties) break;

                // terminate if the next property is not eligible for combining in a multi-pane
                handler = createCurrentPropertyHandler(+1);
                if (!eligableForMultiPane(handler.get())) break;
                _propertyIndex++;
            }
            return multipane;
        }
    case State::SaveHierarchy:
        {
            return new SaveWizardPane(_schema.get(), _root.get(), _filepath, _dirty, this);
        }
    }
    return nullptr;
}

////////////////////////////////////////////////////////////////////

QString WizardEngine::hierarchyPath()
{
    string result;

    if (_state == State::ConstructHierarchy)
    {
        // on the lowest level, show item type and property name
        result = _current->type() + " : " + _schema->properties(_current->type())[_propertyIndex];

        // for higher levels, show only item type
        Item* current = _current->parent();
        while (current)
        {
            result = current->type() + u8" \u2192 " + result;
            current = current->parent();
        }
    }
    return QString::fromStdString(result);
}

////////////////////////////////////////////////////////////////////
