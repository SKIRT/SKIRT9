/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BOOLEANEXPRESSION_HPP
#define BOOLEANEXPRESSION_HPP

#include "Basics.hpp"
#include <functional>

////////////////////////////////////////////////////////////////////

/** The BooleanExpression class allows evaluating a Boolean expression or a conditional value
    expression given the Boolean value for each of the identifiers in the expression.

    An identifier in an expression consists of an arbitrary sequence of (upper- and lowercase)
    letters and digits, starting with a letter. The latter requirement is not verified by the
    current implementation, but is imposed to allow possible future extensions to recognize, for
    example, integer constants. The caller is expected to provide a call-back function to determine
    the Boolean value for each identifier in the expression.

    <B>Boolean expression</B>

    In a Boolean expression, identifiers are combined with the unary negation operator "!", the
    binary OR operator |, the binary AND operator \&, and (balanced) parenthesis to determine the
    ordering of operations. In the absence of parenthesis, the negation operator has the highest
    precedence; otherwise operations are executed from left to right. When mixing OR and AND
    operators it is therefore highly recommended to use parenthesis.

    Boolean expressions are not allowed to contain whitespace. By convention, an empty Boolean
    expression has a value of true.

    <B>Conditional value expression</B>

    The syntax of a conditional value expression can be illustrated as follows:

        <Boolean-expression1>:<value1>; ... ;[<Boolean-expressionN>:]<valueN>

    In other words, there is a sequence of one or more condition-value pairs separated by a
    semicolon, while the condition and the value are separated by a colon. The pairs in the
    conditional value expression are considered from left to right. When the condition for a pair
    evaluates to true (using the Boolean expression rules described above), the corresponding value
    represents the value of the complete expression and further evaluation is terminated. If none
    of the conditions evaluate to true, the value of the expression is the empty string.

    The condition in a pair (and the colon separating it from the value) may be omitted. In this
    case, the condition defaults to true. This makes sense only for the last (or only) pair in the
    sequence, because any remaining pairs will never be considered.

    The value of a conditional value expression is always a string, which may be empty. Conditional
    value expressions are not allowed to contain whitespace except as part of the value. The
    overall expression syntax disallows a value to contain colons and semicolons. This is no
    problem for values that represent a number or a type name, and even for other types of values
    it is most likely not a significant limitation. */
class BooleanExpression final
{
public:
    // ================== Queries ==================

    /** This function evaluates the specified string as a Boolean expression in the format decribed
        in the class header, and returns the result. When evaluating the expression, each
        identifier is replaced by the value returned by the specified callback function. The
        function throws an error if the expression string does not conform to the syntax of a
        Boolean expression. */
    static bool evaluateBoolean(string expression, std::function<bool(string)> isIdentifierTrue);

    /** This function evaluates the specified string as a conditional value expression in the
        format decribed in the class header, and returns the result. When evaluating the
        expression, each identifier is replaced by the value returned by the specified callback
        function. The function throws an error if the expression string does not conform to the
        syntax of a conditional value expression. */
    static string evaluateConditionalValue(string expression, std::function<bool(string)> isIdentifierTrue);
};

////////////////////////////////////////////////////////////////////

#endif
