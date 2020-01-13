/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef XMLREADER_HPP
#define XMLREADER_HPP

#include "Basics.hpp"
#include <fstream>
#include <unordered_map>

////////////////////////////////////////////////////////////////////

/** The XmlReader class implements a parser for the subset of XML that is sufficient to represent
    SMILE schemas and datasets. Specifically, the input stream must be a well-formed XML 1.0
    document with the following additional limitations and caveats.

    - The text must be UTF-8 encoded (upwards compatible with 7-bit ASCII). Byte sequences that
      are not valid UTF-8 result in an error.
    - The presence of a document type definition (DTD) is not supported and results in an error.
    - Only (nested) elements and element attributes are supported. Comments and processing
      instructions (including the XML prolog declaration on the first line) are ignored.
      CDATA sections or text other than white space between elements result in an error.
    - Element and attribute names must start with a 7-bit ASCII letter [a-zA-Z] and can further
      contain letters [a-zA-Z], digits [0-9], and dashes (-). Any other characters are prohibited
      and will result in an error (note that this is a lot more restrictive than standard XML).
    - There is no support for namespaces; a colon in a name results in an error (per the above).
    - Attribute value strings may be surrounded by single or double quotes. To allow attribute
      values to contain both single and double quotes, the single-quote character (') may be
      represented as "&apos;", and the double-quote character (") as "\&quot;".
    - The ampersand character (\&) and the angle brackets (<,>) in attribute values should be
      represented as "&amp;", "&lt;", and "&gt;" respectively. Other than these, entity
      and character references (starting with an \&) are not supported and result in an error.
    - Newline characters are not allowed in attribute values, and a tab is interpreted as a space.

    When an error occurs while opening or parsing the XML data stream, the constructor and methods
    of this class throw a FatalError with an appropriate error message.
*/
class XmlReader final
{
public:
    /** This constructor accepts and retains a reference to the input stream to be parsed and
        initializes the XML parser. It verifies that the stream is ready for reading and consumes
        the UTF-8 BOM if it is present at the beginning of the stream. The second argument provides
        a human readable string to identify the stream in error messages. */
    XmlReader(std::istream& inputStream, string streamName);

    /** This constructor opens the file with the specified path as the input stream to be parsed
        and initializes the XML parser. It verifies that the file is ready for reading and consumes
        the UTF-8 BOM if it is present at the beginning of the file. The file path is used to
        identify the stream in error messages. */
    XmlReader(string filePath);

    /** The copy constructor is deleted because instances of this class should never be copied or
        moved. */
    XmlReader(const XmlReader&) = delete;

    /** The assignment operator is deleted because instances of this class should never be copied
        or moved. */
    XmlReader& operator=(const XmlReader&) = delete;

    // ================== Reading elements ==================

    /** Reads until the next start element within the current element. Returns true when a start
        element was reached, and false when the end element for the current element was reached.
        The current element is the element matching the most recently parsed start element of which
        a matching end element has not yet been reached. When the parser has reached the end
        element, the parent element becomes the current element. */
    bool readNextStartElement();

    /** Reads and ignores any remaining contents of the current element, i.e. any (nested) child
        elements, including the element's end tag, causing the current element's parent element to
        become current. */
    void skipCurrentElement();

    /** Returns the name of the current element, or the empty string if there is none. */
    string elementName() const;

    /** Returns a list with the names of the attributes of the current element. The list is empty
        if there is no current element or if the current element has no attributes. */
    vector<string> attributeNames() const;

    /** Returns the value of the current element's attribute with the specified name. If there is
        no current element or if the current element has no attribute with the given name, an empty
        string is returned. */
    string attributeValue(string name) const;

    // ================== Error handling ==================

    /** Throws an error with the specified error message, augmented with information about the
        current line number in the input stream. This function is used by the XML reader when an
        error occurs, and it may also be used by a client to raise a custom error. */
    void throwError(string message);

    // ================== Private parsing utilities ==================

private:
    /** Peeks at the next character without consuming it. The character is returned as type char
        (which may be signed or unsigned). If there are no more characters in the input stream, the
        function returns a system-dependent value that should differ from regular characters
        (probably 0xFF). */
    char peek();

    /** Gets exactly one character. */
    char get();

    /** Skips exactly one character. */
    void skip();

    /** Skips (a portion of) the UTF-8 byte order marker, if present. */
    void skipBOM();

    /** Skips consecutive white space (space, tab, newline), if any, while counting lines. */
    void skipWhiteSpace();

    /** Skips anything up to and including the given string, while counting lines. The given string
        should not be empty and should not contain newline characters. */
    void skipUpTo(string match);

    /** Gets and returns a name token starting with a letter and further consisting of letters and
        digits. The name token can be terminated by whitespace or by the '=' sign; the terminating
        character is not consumed. If any other characters are encountered, the function throws an
        error. */
    string getName();

    /** Gets and returns an attribute value string surrounded by single or double quotes. Quotes of
        the other type may occur in the string, as well as the escape sequences "&apos;",
        "\&quot;", "&amp;", "&lt;", and "&gt;". If any newline or control characters are
        encountered in the value string, the function throws an error. */
    string getValue();

    // ================== Data members ==================

private:
    std::ifstream _infile;        // the input file, if it is opened in our constructor
    std::istream& _in;            // reference to the input stream (which may be our inout file)
    string _streamName;           // human readable name for use in error messages
    uint64_t _lineNumber{1};      // current line number
    bool _elementIsEmpty{false};  // true if the current element is empty (i.e. there is no separate end-tag)

    struct ElementInfo
    {
        string elementName;                                  // name of the element
        vector<string> attributeNames;                       // attribute names for the element, in order of appearance
        std::unordered_map<string, string> attributeValues;  // attribute key-value pairs for the element
    };
    vector<ElementInfo> _elementStack;  // information for the ancestors of the current element, including itself
};

////////////////////////////////////////////////////////////////////

#endif
