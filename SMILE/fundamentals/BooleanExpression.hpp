/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BOOLEANEXPRESSION_HPP
#define BOOLEANEXPRESSION_HPP

#include "Basics.hpp"
#include <functional>
#include <unordered_set>

////////////////////////////////////////////////////////////////////

/** The BooleanExpression class allows evaluating a Boolean expression with specified values for
    the identifiers in the expression. An identifier in the expression consists of an arbitrary
    sequence of letters and digits. Identifiers are combined with the unary negation operator "!",
    the binary OR operator |, the binary AND operator \&, and (balanced) parenthesis to determine
    the ordering of operations. In the absence of parenthesis, the negation operator has the
    highest precedence; otherwise operations are executed from left to right. When mixing OR and
    AND operators it is therefore highly recommended to use parenthesis. */
class BooleanExpression final
{
public:
    // ================== Queries ==================

    /** This function evaluates the specified string as a Boolean expression in the format decribed
        in the class header, and returns the result. When evaluating the expression, an identifier
        is replaced by the value returned by the specified callback function. The function throws
        an error if the expression string is empty, contains whitespace, or does not conform to the
        syntax of a Boolean expression as described above. */
    static bool evaluate(string expression, std::function<bool(string)> isIdentifierTrue);

    /** This function evaluates the specified string as a Boolean expression in the format decribed
        in the class header, and returns the result. When evaluating the expression, an identifier
        is replaced by true if the identifier is found in the specified set of keywords, and by
        false if it is not. The function throws an error if the expression string is empty,
        contains whitespace, or does not conform to the syntax of a Boolean expression as described
        above. */
    static bool evaluate(string expression, const std::unordered_set<string>& trueIdentifiers);
};

////////////////////////////////////////////////////////////////////

#endif
