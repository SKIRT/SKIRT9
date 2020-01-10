/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LatexHierarchyWriter.hpp"
#include "BoolPropertyHandler.hpp"
#include "DoubleListPropertyHandler.hpp"
#include "DoublePropertyHandler.hpp"
#include "EnumPropertyHandler.hpp"
#include "FatalError.hpp"
#include "IntPropertyHandler.hpp"
#include "Item.hpp"
#include "ItemListPropertyHandler.hpp"
#include "ItemPropertyHandler.hpp"
#include "PropertyHandlerVisitor.hpp"
#include "SchemaDef.hpp"
#include "StringPropertyHandler.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include <unordered_map>

////////////////////////////////////////////////////////////////////

namespace
{
    // This local class offers functions to insert TeX commands into text string of various types.
    // It also holds data stuctures with the replacement pairs; these are initialized "lazily", i.e. only when used.
    class Texify
    {
    private:
        // replacement pairs for characters in regular text (just a small and arbitrary selection)
        const vector<std::pair<string, string>> _replText = {
            // latex escapes
            {"\\", "\\textbackslash{}"},
            {"~", "{\\raise.17ex\\hbox{$\\scriptstyle\\sim$}}"},
            {"#", "\\#"},
            {"$", "\\$"},
            {"%", "\\%"},
            {"&", "\\&"},
            {"_", "\\_"},
            {"{", "\\{"},
            {"}", "\\}"},
            // greek
            {"α", "$\\alpha$"},
            {"β", "$\\beta$"},
            {"γ", "$\\gamma$"},
            {"δ", "$\\delta$"},
            {"ϵ", "$\\epsilon$"},
            {"ε", "$\\varepsilon$"},
            {"ζ", "$\\zeta$"},
            {"η", "$\\eta$"},
            {"θ", "$\\theta$"},
            {"ϑ", "$\\vartheta$"},
            {"κ", "$\\kappa$"},
            {"λ", "$\\lambda$"},
            {"μ", "$\\mu$"},
            {"ν", "$\\nu$"},
            {"ξ", "$\\xi$"},
            {"π", "$\\pi$"},
            {"ρ", "$\\rho$"},
            {"σ", "$\\sigma$"},
            {"τ", "$\\tau$"},
            {"ϕ", "$\\phi$"},
            {"φ", "$\\varphi$"},
            {"χ", "$\\chi$"},
            {"ψ", "$\\psi$"},
            {"ω", "$\\omega$"},
            {"Γ", "$\\Gamma$"},
            {"Δ", "$\\Delta$"},
            {"Θ", "$\\Theta$"},
            {"Λ", "$\\Lambda$"},
            {"Π", "$\\Pi$"},
            {"Σ", "$\\Sigma$"},
            {"Φ", "$\\Phi$"},
            {"Ψ", "$\\Psi$"},
            {"Ω", "$\\Omega$"},
            // accented characters
            {"á", "\\'{a}"},
            {"à", "\\`{a}"},
            {"â", "\\^{a}"},
            {"ä", "\\\"{a}"},
            {"å", "\\r{a}"},
            {"Å", "\\r{A}"},
            {"é", "\\'{e}"},
            {"è", "\\`{e}"},
            {"ê", "\\^{e}"},
            {"ë", "\\\"{e}"},
            {"î", "\\^{\\i}"},
            {"ï", "\\\"{\\i}"},
            {"ó", "\\'{o}"},
            {"ò", "\\`{o}"},
            {"ô", "\\^{o}"},
            {"ö", "\\\"{o}"},
            {"ú", "\\'{u}"},
            {"ù", "\\`{u}"},
            {"û", "\\^{u}"},
            {"ü", "\\\"{u}"},
            {"ç", "\\c{c}"},
            // special symbols
            {"∅", "$\\emptyset$"},
            {"∞", "$\\infty$"},
            {"×", "$\\times$"},
            {"°", "$^{\\circ}$"}};

        // replacement pairs for unit tokens
        const std::unordered_map<string, string> _replUnit = {{"micron", "\\mu\\textrm{m}"},
                                                              {"Angstrom", "\\textrm{\\r{A}}"},
                                                              {"Lsun", "\\textrm{L}_{\\odot}"},
                                                              {"Msun", "\\textrm{M}_{\\odot}"},
                                                              {"deg", "^{\\circ}"},
                                                              {"arcsec", "^{\\prime\\prime}"},
                                                              {"arcsec2", "(^{\\prime\\prime})^{2}"}};

    private:
        // returns pair (number of characters to be removed from start of input string, replacement string)
        // input string should be non-empty!!
        std::pair<size_t, string> nextReplace(string str) const
        {
            for (const auto& find_replace : _replText)
            {
                if (StringUtils::startsWith(str, find_replace.first))
                {
                    return std::make_pair(find_replace.first.size(), find_replace.second);
                }
            }
            return std::make_pair(1, str.substr(0, 1));
        }

    public:
        // texify regular text
        string fromText(string str) const
        {
            // process from left to right to avoid replacing parts of replaced strings
            string result;
            while (!str.empty())
            {
                auto pair = nextReplace(str);
                result += pair.second;
                str.erase(0, pair.first);
            }
            return result;
        }

        // texify regular text after making first character uppercase
        string fromTextUpp(string str) const { return fromText(StringUtils::toUpperFirst(str)); }

        // texify string representing a floating point value in SMILE format
        string fromDouble(string str) const
        {
            // get the segments
            vector<string> segments = StringUtils::split(StringUtils::squeeze(str), " ");
            if (segments.empty() || segments.size() > 2) return str;

            // treat the number (first handle the exponent since the special symbol replacements also contain 'e')
            string result = segments[0];
            if (StringUtils::contains(result, "e"))
            {
                result = StringUtils::replace(result, "e", "\\times 10^{");
                result += "}";
            }
            result = StringUtils::replace(result, "∅", "\\emptyset");
            result = StringUtils::replace(result, "∞", "\\infty");
            result = StringUtils::replace(result, "-", "\\textrm{-}");

            // treat the unit if present
            if (segments.size() > 1)
            {
                // split tokens separated by slashes
                vector<string> tokens = StringUtils::split(segments[1], "/");

                for (string& token : tokens)  // writable reference
                {
                    // replace the complete token if it is in the list of specials
                    // otherwise, put in in upright type (except for an exponent at the end)
                    if (_replUnit.count(token))
                        token = _replUnit.at(token);
                    else
                    {
                        auto lastChar = token.back();
                        if (lastChar == '2' || lastChar == '3' || lastChar == '4')
                            token = "\\textrm{" + token.substr(token.size() - 1) + "}^{" + lastChar + "}";
                        else
                            token = "\\textrm{" + token + "}";
                    }
                }
                result += "\\:" + StringUtils::join(tokens, "/");
            }

            // return result as math expression
            return "$" + result + "$";
        }

        // texify string representing a list of floating point values in SMILE format
        string fromDoubleList(string str) const
        {
            vector<string> result;
            for (string item : StringUtils::split(str, ",")) result.push_back(fromDouble(item));
            return StringUtils::join(result, ", ");
        }
    };
}

////////////////////////////////////////////////////////////////////

namespace
{
    class LatexWriter : PropertyHandlerVisitor
    {
    private:
        const SchemaDef* _schema;
        std::ofstream _out;
        size_t _indent{0};
        const Texify tex{};

    public:
        LatexWriter(Item* item, const SchemaDef* schema, string filePath, string dataset, string producer)
            // remember schema and open output stream
            : _schema(schema), _out(System::ofstream(filePath))
        {
            // write document preamble
            _out << "\\documentclass[10pt,english]{article}\n";
            _out << "\\usepackage[landscape,a4paper]{geometry}\n";
            _out << "\\geometry{verbose,tmargin=2cm,bmargin=2cm,lmargin=3cm,rmargin=2cm}\n";
            _out << "\\pagestyle{empty}\n";
            _out << "\\setlength{\\parskip}{\\bigskipamount}\n";
            _out << "\\setlength{\\parindent}{0pt}\n";
            _out << "\\usepackage{babel}\n";
            _out << "\\begin{document}\n";
            _out << "\n";

            // write header text
            _out << "\\section*{" << tex.fromText(dataset) << " (" << tex.fromText(schema->schemaTitle()) << ")}\n";
            _out << "Generated by " << tex.fromText(!producer.empty() ? producer : "SMILE") << " on "
                 << System::timestamp() << "\\\\\n";
            _out << "\\copyright Astronomical Observatory, Ghent University\\\\\n";
            _out << "\n";

            // write the first label/value line for the top-level item
            _out << tex.fromTextUpp(schema->title(item->type())) << "\\\\\n";

            // recursively write all properties of the top-level item and its children
            writeProperties(item);

            // write document end
            _out << "\n";
            _out << "\\end{document}\n";
        }

        void writeProperties(Item* item)
        {
            _indent++;

            // handle all properties of the item
            for (const string& property : _schema->properties(item->type()))
            {
                auto handler = _schema->createPropertyHandler(item, property, nullptr);
                // distribute to write methods depending on property type (using visitor pattern)
                handler->acceptVisitor(this);
            }

            _indent--;
        }

        void visitPropertyHandler(StringPropertyHandler* handler) override
        {
            // write a label/value line for this property
            writeIndent();
            _out << tex.fromTextUpp(handler->title()) << ": " << tex.fromText(handler->value()) << "\\\\\n";
        }

        void visitPropertyHandler(BoolPropertyHandler* handler) override
        {
            // write a label/value line for this property
            writeIndent();
            _out << tex.fromTextUpp(handler->title()) << ": " << (handler->value() ? "yes" : "no") << "\\\\\n";
        }

        void visitPropertyHandler(IntPropertyHandler* handler) override
        {
            // write a label/value line for this property
            writeIndent();
            _out << tex.fromTextUpp(handler->title()) << ": $" << StringUtils::toString(handler->value()) << "$\\\\\n";
        }

        void visitPropertyHandler(EnumPropertyHandler* handler) override
        {
            // write a label/value line for this property
            writeIndent();
            _out << tex.fromTextUpp(handler->title()) << ": " << tex.fromText(handler->titleForValue()) << "\\\\\n";
        }

        void visitPropertyHandler(DoublePropertyHandler* handler) override
        {
            // write a label/value line for this property
            writeIndent();
            _out << tex.fromTextUpp(handler->title()) << ": " << tex.fromDouble(handler->toString(handler->value()))
                 << "\\\\\n";
        }

        void visitPropertyHandler(DoubleListPropertyHandler* handler) override
        {
            // write a label/value line for this property
            writeIndent();
            _out << tex.fromTextUpp(handler->title()) << ": " << tex.fromDoubleList(handler->toString(handler->value()))
                 << "\\\\\n";
        }

        void visitPropertyHandler(ItemPropertyHandler* handler) override
        {
            if (handler->value())
            {
                // write a label/value line for this property
                writeIndent();
                _out << tex.fromTextUpp(handler->title()) << ": "
                     << tex.fromText(_schema->title(handler->value()->type())) << "\\\\\n";

                // handle the properties of the item pointed to by the property
                writeProperties(handler->value());
            }
        }

        void visitPropertyHandler(ItemListPropertyHandler* handler) override
        {
            // loop over all items in the list held by this property
            size_t index = 0;
            for (Item* item : handler->value())
            {
                index++;

                // write a label/value line for this item in the list
                writeIndent();
                _out << "Item \\#" << std::to_string(index) << " in " << tex.fromText(handler->title())
                     << " list: " << tex.fromText(_schema->title(item->type())) << "\\\\\n";

                // handle the properties of the item pointed to by this item in the list
                writeProperties(item);
            }
        }

        void writeIndent()
        {
            for (size_t i = 0; i != _indent; ++i) _out << ".\\quad{}";
        }
    };
}

////////////////////////////////////////////////////////////////////

void LatexHierarchyWriter::write(Item* item, const SchemaDef* schema, string filePath, string dataset, string producer)
{
    LatexWriter writer(item, schema, filePath, dataset, producer);
}

////////////////////////////////////////////////////////////////////
