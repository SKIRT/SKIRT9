/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef XMLWRITER_HPP
#define XMLWRITER_HPP

#include "Basics.hpp"
#include <fstream>

////////////////////////////////////////////////////////////////////

/** The XmlWriter class facilates writing the subset of XML that is sufficient to represent
    SMILE schemas and datasets. The output stream is garantueed to be a well-formed XML 1.0
    document that also conforms to the limitations described for the XmlReader class, so that
    it can be parsed by that class.

    When an error occurs while writing the XML data stream, the constructor and methods
    of this class throw a FatalError with an appropriate error message.
*/
class XmlWriter final
{
public:
    /** This constructor accepts and retains a reference to the output stream to be written, and
        initializes the XML writer. The second argument provides a human readable string to
        identify the stream in error messages. */
    XmlWriter(std::ostream& outputStream, string streamName);

    /** This constructor creates and opens the file with the specified path as the output stream to
        be written, and initializes the XML writer. The file path is used to identify the stream in
        error messages. */
    XmlWriter(string filePath);

    /** The copy constructor is deleted because instances of this class should never be copied or
        moved. */
    XmlWriter(const XmlWriter&) = delete;

    /** The assignment operator is deleted because instances of this class should never be copied
        or moved. */
    XmlWriter& operator=(const XmlWriter&) = delete;

    // ================== Public facilities  ==================

    /** Writes a document start with the appropriate XML version and encoding information. This
        function must be called once before any of the other functions are called. */
    void writeStartDocument();

    /** Writes the specified text as XML comment. The text must not contain the forbidden sequence
        "--" or end with "-". Note that XML does not provide any way to escape "-" in a comment. */
    void writeComment(string text);

    /** Writes a start element with the specified name. Subsequent calls to writeAttribute() will
        add attributes to this element. */
    void writeStartElement(string name);

    /** Writes an attribute with the specified name and value. This function can only be called
        after writeStartElement() before any further content is written. */
    void writeAttribute(string name, string value);

    /** Closes the previous matching start element. For a particular XML document, the number of
        invocations of this function must exactly match the number of invocations of the
        writeStartElement() function. This is verified in the writeEndDocument() function. */
    void writeEndElement();

    /** Verifies that the number of writeStartElement() and writeEndElement() invocations match,
        and that the output stream is in a non-error state. If not, a FatalError is thrown. This
        function must be called once after any of the other functions are called. */
    void writeEndDocument();

    // ================== Private utilities ==================

private:
    /** Writes \f$4\times N\f$ spaces where \f$N\f$ is the current element nesting depth. */
    void writeIndentation();

    // ================== Data members ==================

private:
    std::ofstream _outfile;        // the output file, if it is opened in our constructor
    std::ostream& _out;            // reference to the output stream (which may be our output file)
    string _streamName;            // human readable name for use in error messages
    vector<string> _elementNames;  // names of the ancestors of the current element, including itself
    bool _starting{false};         // becomes true while a start element tag is being written and needs closing
};

////////////////////////////////////////////////////////////////////

#endif
