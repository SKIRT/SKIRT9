/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "XmlReader.hpp"
#include "FatalError.hpp"
#include "StringUtils.hpp"
#include "System.hpp"

////////////////////////////////////////////////////////////////////

XmlReader::XmlReader(std::istream& inputStream, string streamName) : _in(inputStream), _streamName(streamName)
{
    // Verify that the stream is operational
    if (_in.peek() == EOF) throw FATALERROR("Can't read from XML input stream or stream is empty: " + _streamName);

    // Skip the UTF-8 BOM, if present
    skipBOM();
}

////////////////////////////////////////////////////////////////////

XmlReader::XmlReader(string filePath) : _in(_infile), _streamName(filePath)
{
    // Open input file and verify that it is operational
    _infile = System::ifstream(filePath);
    if (_in.peek() == EOF) throw FATALERROR("Can't open XML input file or file is empty: " + _streamName);

    // Skip the UTF-8 BOM, if present
    skipBOM();
}

////////////////////////////////////////////////////////////////////

bool XmlReader::readNextStartElement()
{
    // If the current element is empty, we need not look for an end-tag
    if (_elementIsEmpty)
    {
        _elementIsEmpty = false;
        _elementStack.pop_back();
        return false;  // >>>>>>>>>>>>>>>>>>> exit point <<<<<<<<<<<<<<<<<
    }

    // Loop to allow skipping comments and processing instructions
    while (true)
    {
        // Since we don't allow text contents, the first non-whitespace character must be "<"
        skipWhiteSpace();
        if (get() != '<')
        {
            throwError("Encountered unsupported text contents while looking for next element");
        }

        // The character(s) after the '<" determine what we're looking at
        char c = peek();
        if (c == '?')  // *** processing instruction -> ignore
        {
            skip();
            skipUpTo("?>");
        }
        else if (c == '!')
        {
            skip();
            if (get() == '-' && get() == '-')  // *** comments -> ignore
            {
                skipUpTo("-->");
            }
            else  // *** something unsupported or plain wrong
            {
                throwError("Encountered unsupported construct (DOCTYPE or CDATA?) while looking for next element");
            }
        }
        else if (c == '/')  // *** end tag
        {
            // Parse element name and closing bracket
            skip();
            string elemName = getName();
            skipWhiteSpace();
            if (get() != '>')
            {
                throwError("Encountered invalid element end-tag");
            }

            // Verify element nesting and pop the current element from the stack
            if (elementName() != elemName)
            {
                throwError("Element end tag '" + elemName + "' does not match start tag '" + elementName() + "'");
            }
            _elementStack.pop_back();
            return false;  // >>>>>>>>>>>>>>>>>>> exit point <<<<<<<<<<<<<<<<<
        }
        else  // *** start tag, empty-element tag, or error (caught later)
        {
            // Parse and store the element name
            ElementInfo info;
            info.elementName = getName();
            _elementIsEmpty = false;

            // Parse attribute key-value pairs and closing bracket (or slash-bracket)
            while (true)
            {
                // Check for closing bracket, possibly preceded by slash
                skipWhiteSpace();
                if (peek() == '/')
                {
                    skip();
                    _elementIsEmpty = true;
                }
                if (peek() == '>')
                {
                    skip();
                    break;
                }

                // Parse and store attribute key and value
                string name = getName();
                skipWhiteSpace();
                if (get() != '=') throwError("Expected equal sign after attribute name");
                skipWhiteSpace();
                string value = getValue();
                info.attributeNames.push_back(name);
                info.attributeValues.emplace(name, value);
            }

            // Push this new element onto the stack
            _elementStack.emplace_back(std::move(info));
            return true;  // >>>>>>>>>>>>>>>>>>> exit point <<<<<<<<<<<<<<<<<
        }
    }

    // This point is never reached; the loop always exits with a return or a throw, but the compiler can't know this
    return false;
}

////////////////////////////////////////////////////////////////////

void XmlReader::skipCurrentElement()
{
    size_t depth = 1;
    while (depth)
    {
        if (readNextStartElement())
            ++depth;
        else
            --depth;
    }
}

////////////////////////////////////////////////////////////////////

string XmlReader::elementName() const
{
    if (_elementStack.empty()) return string();
    return _elementStack.back().elementName;
}

////////////////////////////////////////////////////////////////////

vector<string> XmlReader::attributeNames() const
{
    if (_elementStack.empty()) return vector<string>();
    return _elementStack.back().attributeNames;
}

////////////////////////////////////////////////////////////////////

string XmlReader::attributeValue(string name) const
{
    if (!_elementStack.empty())
    {
        auto& map = _elementStack.back().attributeValues;
        auto pair = map.find(name);
        if (pair != map.cend()) return pair->second;
    }
    return string();
}

////////////////////////////////////////////////////////////////////

void XmlReader::throwError(string message)
{
    throw FATALERROR("Error in XML stream '" + _streamName + "' at line " + std::to_string(_lineNumber) + "\n"
                     + message);
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

char XmlReader::peek()
{
    return static_cast<char>(_in.peek());
}

////////////////////////////////////////////////////////////////////

char XmlReader::get()
{
    int c = _in.get();
    if (c == EOF) throwError("Reached end of input stream unexpectedly");
    return static_cast<char>(c);
}

////////////////////////////////////////////////////////////////////

void XmlReader::skip()
{
    int c = _in.get();
    if (c == EOF) throwError("Reached end of input stream unexpectedly");
}

////////////////////////////////////////////////////////////////////

void XmlReader::skipBOM()
{
    const char utf8_BOM1 = static_cast<char>(0xEFu);
    const char utf8_BOM2 = static_cast<char>(0xBBu);
    const char utf8_BOM3 = static_cast<char>(0xBFu);

    if (peek() != utf8_BOM1) return;
    skip();
    if (peek() != utf8_BOM2) return;
    skip();
    if (peek() != utf8_BOM3) return;
    skip();
}

////////////////////////////////////////////////////////////////////

void XmlReader::skipWhiteSpace()
{
    // We assume that when entering this function, the stream is not positioned
    // between a carriage return and a line feed
    bool previousCR = false;

    while (true)
    {
        switch (peek())
        {
            case ' ':
            case '\t': previousCR = false; break;

            case 0x0D:  // carriage return
                _lineNumber++;
                previousCR = true;
                break;

            case 0x0A:  // line feed
                if (!previousCR) _lineNumber++;
                previousCR = false;
                break;

            default: return;
        }
        skip();
    }
}

////////////////////////////////////////////////////////////////////

void XmlReader::skipUpTo(string match)
{
    // We assume that when entering this function, the stream is not positioned
    // between a carriage return and a line feed
    bool previousCR = false;

    // The number of characters already matched; we assume match does not contain CR or LF
    size_t matched = 0;

    while (true)
    {
        char c = get();
        if (c == 0x0D)  // carriage return
        {
            _lineNumber++;
            previousCR = true;
            matched = 0;
        }
        else if (c == 0x0A)  // line feed
        {
            if (!previousCR) _lineNumber++;
            previousCR = false;
            matched = 0;
        }
        else if (c == match[matched])
        {
            matched++;
            if (matched == match.length()) return;
            previousCR = false;
        }
        else
        {
            previousCR = false;
            matched = 0;
        }
    }
}

////////////////////////////////////////////////////////////////////

namespace
{
    bool isLetter(char c) { return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'); }

    bool isLetterOrDigitOrDash(char c)
    {
        return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '-';
    }

    bool isNameTerminator(char c)
    {
        return c == ' ' || c == '\t' || c == 0x0D || c == 0x0A || c == '=' || c == '>' || c == '/';
    }

    bool isControlCharacter(char c)  // including newline characters
    {
        return (c >= 0x00 && c < ' ' && c != '\t') || c == 0x7F;
    }
}

////////////////////////////////////////////////////////////////////

string XmlReader::getName()
{
    // Get the first character, which must be a letter
    char c = get();
    if (!isLetter(c))
    {
        throwError("Element or attribute name must start with a letter");
    }
    string name(1, c);

    // Get subsequent characters, which may be letters or digits
    while (true)
    {
        char c = peek();
        if (isNameTerminator(c)) return name;
        if (!isLetterOrDigitOrDash(c))
        {
            throwError("Element or attribute name must consist of letters and digits");
        }
        name += c;
        skip();
    }
}

////////////////////////////////////////////////////////////////////

string XmlReader::getValue()
{
    // Get the string delimiter
    char quote = get();
    if (quote != '\'' && quote != '\"') throwError("Expected single or double quote at start of attribute value");

    // Get the raw contents of the string
    string value;
    while (true)
    {
        char c = get();
        if (c == quote)  // check for terminator
        {
            break;
        }
        if (c == '\t')  // replace tab by space
        {
            c = ' ';
        }
        if (isControlCharacter(c))  // disallow control characters
        {
            throwError("Newline and control characters are not allowed within an attribute value");
        }
        value += c;
    }

    // Substitute the escape sequences
    value = StringUtils::replace(value, "&apos;", "\'");
    value = StringUtils::replace(value, "&quot;", "\"");
    value = StringUtils::replace(value, "&lt;", "<");
    value = StringUtils::replace(value, "&gt;", ">");
    value = StringUtils::replace(value, "&amp;", "&");  // do this one last to avoid introducing new escape sequences

    return value;
}

////////////////////////////////////////////////////////////////////
