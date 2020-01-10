/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "XmlWriter.hpp"
#include "FatalError.hpp"
#include "StringUtils.hpp"
#include "System.hpp"

////////////////////////////////////////////////////////////////////

XmlWriter::XmlWriter(std::ostream& outputStream, string streamName) : _out(outputStream), _streamName(streamName)
{
    // Verify that the stream is operational
    if (!_out) throw FATALERROR("Can't write to XML output stream: " + _streamName);
}

////////////////////////////////////////////////////////////////////

XmlWriter::XmlWriter(string filePath) : _out(_outfile), _streamName(filePath)
{
    // Open output file and verify that it is operational
    _outfile = System::ofstream(filePath);
    if (!_out) throw FATALERROR("Can't create XML output file: " + _streamName);
}

////////////////////////////////////////////////////////////////////

namespace
{
    bool isLetter(char c) { return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'); }

    bool isLetterOrDigitOrDash(char c)
    {
        return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '-';
    }

    bool isControlCharacter(char c)  // including newline characters
    {
        return (c >= 0x00 && c < ' ' && c != '\t') || c == 0x7F;
    }

    // a valid element/attribute name starts with a letter and contains letters, digits and dashes
    bool isValidName(string name)
    {
        return !name.empty() && isLetter(name[0]) && std::all_of(name.cbegin() + 1, name.cend(), isLetterOrDigitOrDash);
    }

    // a valid attribute value does not contain control characters
    bool isValidValue(string value) { return std::none_of(value.cbegin(), value.cend(), isControlCharacter); }

    // a valid comment does not contain control characters or the sequence "--", and does not end in "-"
    bool isValidComment(string text)
    {
        return isValidValue(text) && text.find("--") == string::npos && (text.empty() || text.back() != '-');
    }

    string escapeSpecialChars(string text)
    {
        text = StringUtils::replace(text, "&", "&amp;");  // do this one first to avoid escaping escape sequences
        text = StringUtils::replace(text, "\'", "&apos;");
        text = StringUtils::replace(text, "\"", "&quot;");
        text = StringUtils::replace(text, "<", "&lt;");
        text = StringUtils::replace(text, ">", "&gt;");
        return text;
    }
}

////////////////////////////////////////////////////////////////////

void XmlWriter::writeStartDocument()
{
    if (!_elementNames.empty())
        throw FATALERROR("Document start should be written before other content for XML output stream: " + _streamName);

    _out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
}

////////////////////////////////////////////////////////////////////

void XmlWriter::writeComment(string text)
{
    // if there is an open start element tag, close it
    if (_starting)
    {
        _out << ">\n";
        _starting = false;
    }

    if (!isValidComment(text)) throw FATALERROR("Invalid comment text for XML output stream: " + _streamName);

    // write comment line with indentation prefix
    writeIndentation();
    _out << "<!--" << escapeSpecialChars(text) << "-->\n";
}

////////////////////////////////////////////////////////////////////

void XmlWriter::writeStartElement(string name)
{
    // if there is an open start element tag, close it
    if (_starting)
    {
        _out << ">\n";
    }

    if (!isValidName(name))
        throw FATALERROR("Invalid element name '" + name + "' for XML output stream: " + _streamName);

    // start a new element with indentation prefix
    writeIndentation();
    _out << "<" << name;
    _starting = true;

    // increase the nesting level
    _elementNames.push_back(name);
}

////////////////////////////////////////////////////////////////////

void XmlWriter::writeAttribute(string name, string value)
{
    if (!_starting)
        throw FATALERROR("Can't write attribute when no element tag is open for XML output stream: " + _streamName);

    if (!isValidName(name))
        throw FATALERROR("Invalid attribute name '" + name + "' for XML output stream: " + _streamName);
    if (!isValidValue(value))
        throw FATALERROR("Invalid value for attribute '" + name + "' for XML output stream: " + _streamName);

    // write attribute segment
    _out << " " << name << "=\"" << escapeSpecialChars(value) << "\"";
}

////////////////////////////////////////////////////////////////////

void XmlWriter::writeEndElement()
{
    if (_elementNames.empty())
        throw FATALERROR("End element has no matching start element for XML output stream: " + _streamName);

    // remember the element name and decrease the nesting level
    string name = _elementNames.back();
    _elementNames.pop_back();

    // if there is an open start element tag, close it as an empty element
    if (_starting)
    {
        _out << "/>\n";
        _starting = false;
    }
    else
    {
        // write end element tag with indentation prefix
        writeIndentation();
        _out << "</" << name << ">\n";
    }
}

////////////////////////////////////////////////////////////////////

void XmlWriter::writeEndDocument()
{
    if (!_out) throw FATALERROR("An error ocurred while writing to XML output stream: " + _streamName);
    if (!_elementNames.empty())
        throw FATALERROR("The number of start and end elements does not match for XML output stream: " + _streamName);
}

////////////////////////////////////////////////////////////////////

void XmlWriter::writeIndentation()
{
    auto count = _elementNames.size();
    while (count--) _out << "    ";
}

////////////////////////////////////////////////////////////////////
