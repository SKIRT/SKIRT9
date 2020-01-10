/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BooleanExpression.hpp"
#include "FatalError.hpp"
#include "StringUtils.hpp"
#include <sstream>

////////////////////////////////////////////////////////////////////

namespace
{
    bool isLetterOrDigit(int c) { return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9'); }

    class BooleanExpressionParser final
    {
    private:
        std::stringstream _in;
        std::function<bool(string)> _isTrue;

    public:
        BooleanExpressionParser(string expression, std::function<bool(string)> isIdentifierTrue)
            : _in(expression), _isTrue(isIdentifierTrue)
        {}

        bool evaluate()
        {
            bool result = expression();
            if (_in.get() != EOF) throw FATALERROR("Invalid Boolean expression");
            return result;
        }

    private:
        bool expression()  // T, T*T, T*T*T, ...   where * is & or |
        {
            // process the first term
            bool result = term();

            while (_in.peek() == '&' || _in.peek() == '|')
            {
                // get the operator and process the subsequent term accordingly
                if (_in.get() == '&')
                    result &= term();
                else
                    result |= term();
            }
            return result;
        }

        bool term()  // (E), !T, identifier
        {
            bool result = false;

            // get the next character
            int c = _in.get();
            if (c == EOF) throw FATALERROR("Invalid Boolean expression");

            // parenthesis
            if (c == '(')
            {
                result = expression();
                if (_in.get() != ')') throw FATALERROR("Invalid Boolean expression");
            }

            // negation
            else if (c == '!')
            {
                result = !term();
            }

            // identifier
            else if (isLetterOrDigit(c))
            {
                // gather the complete identifier without consuming the next character
                string identifier(1, c);
                while (isLetterOrDigit(_in.peek())) identifier += _in.get();

                // get its value from the call-back function
                result = _isTrue(identifier);
            }

            // unsupported character
            else
                throw FATALERROR("Invalid Boolean expression");

            return result;
        }
    };
}

////////////////////////////////////////////////////////////////////

bool BooleanExpression::evaluateBoolean(string expression, std::function<bool(string)> isIdentifierTrue)
{
    if (expression.empty()) return true;
    BooleanExpressionParser parser(expression, isIdentifierTrue);
    return parser.evaluate();
}

////////////////////////////////////////////////////////////////////

string BooleanExpression::evaluateConditionalValue(string expression, std::function<bool(string)> isIdentifierTrue)
{
    // loop over all pairs in the expression
    for (string pair : StringUtils::split(expression, ";"))
    {
        auto splitpair = StringUtils::split(pair, ":");
        if (splitpair.size() > 2) throw FATALERROR("Invalid conditional value expression");

        // if there is no colon, this pair becomes the result
        if (splitpair.size() == 1) return pair;

        // if there is a colon, evaluate the condition
        if (evaluateBoolean(splitpair[0], isIdentifierTrue)) return splitpair[1];
    }
    return string();
}

////////////////////////////////////////////////////////////////////
