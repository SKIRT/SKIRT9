/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Block.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

// lines will be wrapped before reaching the margin
const static size_t MARGIN = 100;

////////////////////////////////////////////////////////////////////

Block::Block(const vector<string>& block)
{
    _block = block;
}

////////////////////////////////////////////////////////////////////

Block::Block(const vector<string>& chunk, size_t first, size_t last)
{
    for (size_t index = first; index <= last; index++)
    {
        _block.push_back(chunk[index]);
    }
}

////////////////////////////////////////////////////////////////////

// this function assumes that _block contains at least one line,
// that the first line starts with /** and that the last line ends with */
vector<string> Block::streamlined()
{
    vector<string> result;

    // get the portion of the first line before the slash
    string prefix = _block[0].substr(0, _block[0].find('/'));

    // initialize the first output line including the leading /**
    string outputline = prefix + "/**";

    // indicate that we can't output an empty line at the start of the block
    bool canWriteEmpty = false;

    // process each input line in turn
    for (size_t index = 0; index < _block.size(); index++)
    {
        // replace consecutive whitespace by a single space
        string inputline = StringUtils::squeeze(_block[index]);

        // remove the leading /** on the first line and the optional * on subsequent lines (including the space)
        if (index == 0) inputline.erase(0, 2);
        if (StringUtils::startsWith(inputline, "* ")) inputline.erase(0, 2);
        if (StringUtils::startsWith(inputline, "*") && !StringUtils::startsWith(inputline, "*/")) inputline.erase(0, 1);

        // if needed, insert a space before the trailing */
        if (index == _block.size() - 1 && inputline.size() > 2 && inputline[inputline.size() - 3] != ' ')
            inputline.insert(inputline.size() - 2, 1, ' ');

        // the first empty input line in a sequence triggers a single empty output line
        if (inputline.empty())
        {
            if (canWriteEmpty)
            {
                if (!StringUtils::squeeze(outputline).empty()) result.push_back(outputline);
                result.push_back("");
                outputline = prefix + "   ";
                canWriteEmpty = false;
            }
        }
        else
        {
            canWriteEmpty = true;

            // split the input line into words and process each word in turn
            for (const string& word : StringUtils::split(inputline, " "))
            {
                // if the output line is full (with one space to spare), flush the line and start a new one
                if (outputline.size() + word.size() > MARGIN - 2)
                {
                    result.push_back(outputline);
                    outputline = prefix + "   ";
                }
                // add the word to the output line
                outputline += " " + word;
            }
        }
    }

    // flush the last output line if nonempty
    if (!StringUtils::squeeze(outputline).empty()) result.push_back(outputline);
    return result;
}

////////////////////////////////////////////////////////////////////
