/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NAMEMANAGER_HPP
#define NAMEMANAGER_HPP

#include "Basics.hpp"
#include <stack>
#include <unordered_set>

////////////////////////////////////////////////////////////////////

/** The NameManager class assists an engine constructing or outputting a SMILE dataset with
    implementing conditional options and values. The overall goal is to dynamically adjust
    displayed options and default values to choices made earlier in the configuration process.
    These earlier choices may represent various aspects of the configuration, including the level
    of expertise selected by the user, some overall mode configured early on, or the value of the
    option just preceding this one. All the metadata needed to control these capabilities are
    defined as part of the SMILE schema.

    Through a NameManager object, an engine processing a SMILE dataset maintains two sets of names
    with a different scope. The global set contains names that stay around for the duration of the
    configuration session. The local set contains names that stay around only while considering the
    properties of a particular type; a fresh local set is created for each type instance. By
    definition, the global set contains names that start with an uppercase letter, and the local
    set contains names that start with a lowercase letter. As a result, the two sets form disjoint
    namespaces. After the first character (which must be a letter), a name can include any mix of
    letters and digits (and no other characters).

    Names are never removed from a set (although the current local set is regularly discarded and
    replaced by a new one as a whole). Names are inserted in each of these sets as the dataset gets
    processed according to the rules set forth by the SMILE specification. Boolean expressions and
    conditional value expressions can then be evaluated as they occur in type and property
    attributes, replacing names contained in one of the sets by true, and other names by false.

    To deal with the nested nature of a SMILE dataset, this class actually maintains a stack of
    local name sets. A fresh local name set can be be pushed onto the stack just before the engine
    starts handling a compound property, and popped from the stack when the engine has completed
    handling the property. */
class NameManager
{
    // ================== Constructing and managing name sets as a whole ==================

public:
    /** The constructor creates a new name manager with both global and local name sets in a fresh
        state. The global name set contains "True" and the local name set contains "true". */
    NameManager();

    /** This functions clears the local name set stack and initializes both global and local name
        sets to their fresh state. After this function completes, the global name set contains
        "True" and the local name set contains "true". */
    void clearAll();

    /** This functions pushes a fresh local name set onto the stack, hiding the previous local name
        set. After this function completes, the local name set contains "true". The global name set
        is not affected. */
    void pushLocal();

    /** This functions pops the local name set from the stack, discarding it, and unhiding the
        local name set that was most recently pushed. The global name set is not affected. A
        popLocal() invocation that does not match a previous pushLocal() invocation causes
        undefined behavior. */
    void popLocal();

    // ================== Inserting names ==================

public:
    /** Adds the specified name to the global or local name set depending on the capitalization of
        the first character. If the first character is an uppercase letter, the name is inserted
        into the global name set. If the first character is a lowercase letter, the name is
        inserted into the local name set. If the first character is not a letter, or if the string
        contains characters other than letters and digits, a fatal error is thrown. */
    void insert(string name);

    /** Adds the specified names to the appropriate name set(s) as described for the insert(string)
        function. */
    void insert(const vector<string>& names);

    /** Adds the names provided in the specified conditional value expression to the appropriate
        name set(s). The conditional value expression is evaluated as described for the
        evaluateConditionalValue() function. The result is interpreted as a comma-separated list of
        names. Each of these names is inserted as described for the insert() function. */
    void insertFromConditionalValue(string nameExpression);

    /** Adds the names provided in the specified conditional value expressions to the appropriate
        name set(s) as described for the insertFromConditionalValue(string) function. */
    void insertFromConditionalValue(const vector<string>& nameExpressions);

    // ================== Evaluating expressions ==================

    /** This function evaluates the specified string as a Boolean expression in the format decribed
        in the class header, and returns the result. When evaluating the expression, each
        identifier is replaced by true if the corresponding name is in the global or the local set,
        and by false if it is not. The function throws an error if the expression string does not
        conform to the syntax of a Boolean expression. */
    bool evaluateBoolean(string expression) const;

    /** This function evaluates the specified string as a conditional value expression in the
        format decribed in the class header, and returns the result. When evaluating the
        expression, each identifier is replaced by true if the corresponding name is in the global
        or the local set, and by false if it is not. The function throws an error if the expression
        string does not conform to the syntax of a conditional value expression. */
    string evaluateConditionalValue(string expression) const;

    // ================== Data members ==================

private:
    std::unordered_set<string> _globalSet;
    std::stack<std::unordered_set<string>> _localSetStack;
};

////////////////////////////////////////////////////////////////////

#endif
