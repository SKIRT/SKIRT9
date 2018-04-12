/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BooleanExpression.hpp"
#include "FatalError.hpp"
#include <sstream>

////////////////////////////////////////////////////////////////////

namespace
{
    bool isLetterOrDigit(int c)
    {
        return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z')  || (c >= '0' && c <= '9');
    }

    class BooleanExpressionParser final
    {
    private:
        std::stringstream _in;
        std::function<bool(string)> _isTrue;

    public:
        BooleanExpressionParser(string expression, std::function<bool(string)> isIdentifierTrue)
            :_in(expression), _isTrue(isIdentifierTrue) { }

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
                if (_in.get() == '&') result &= term();
                else result |= term();
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

                // check presence in the identifier list
                result = _isTrue(identifier);
            }

            // unsupported character
            else throw FATALERROR("Invalid Boolean expression");;

            return result;
        }
    };
}

////////////////////////////////////////////////////////////////////

bool BooleanExpression::evaluate(string expression, std::function<bool(string)> isIdentifierTrue)
{
    BooleanExpressionParser parser(expression, isIdentifierTrue);
    return parser.evaluate();
}

////////////////////////////////////////////////////////////////////

bool BooleanExpression::evaluate(string expression, const std::unordered_set<string>& trueIdentifiers)
{
    BooleanExpressionParser parser(expression, [&trueIdentifiers] (string identifier)
    {
        return trueIdentifiers.count(identifier) > 0;
    });
    return parser.evaluate();
}

////////////////////////////////////////////////////////////////////
