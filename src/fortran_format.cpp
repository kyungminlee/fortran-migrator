#include "fortran_format.h"

#include <algorithm>
#include <cctype>

namespace fmig {

FortranForm detect_form(std::string_view extension) {
    std::string ext(extension);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    if (ext == ".f90" || ext == ".f95" || ext == ".f03" || ext == ".f08")
        return FortranForm::Free;
    // .f, .for, .ftn, and uppercase variants (.F, .F90, etc.)
    // .F90 etc. are free-form with preprocessing; detect by lowercase form.
    return FortranForm::Fixed;
}

bool is_fixed_form_comment(std::string_view line) {
    if (line.empty())
        return true; // Blank lines are treated as comments in fixed-form.
    char c = line[0];
    return c == 'C' || c == 'c' || c == '*' || c == '!';
}

bool is_fixed_form_continuation(std::string_view line) {
    // A continuation line has a non-blank, non-zero character in column 6
    // (0-indexed column 5), and columns 1-5 are blank or numeric label.
    if (line.size() < 6)
        return false;
    if (is_fixed_form_comment(line))
        return false;
    char col6 = line[5];
    return col6 != ' ' && col6 != '0';
}

std::string_view fixed_form_statement(std::string_view line) {
    if (line.size() <= 6)
        return {};
    size_t end = std::min(line.size(), size_t(72));
    return line.substr(6, end - 6);
}

std::string reformat_fixed_line(std::string_view line,
                                char continuation_char) {
    // If line fits in 72 columns, return as-is.
    if (line.size() <= 72)
        return std::string(line);

    // Split the statement portion across continuation lines.
    // Columns 1-6 of the first line stay as-is.
    std::string result;
    std::string_view prefix = line.substr(0, 6);
    std::string_view body = line.substr(6);

    constexpr size_t stmt_width = 66; // columns 7-72
    bool first = true;
    while (!body.empty()) {
        size_t chunk = std::min(body.size(), stmt_width);
        if (first) {
            result += prefix;
            first = false;
        } else {
            result += "\n     ";
            result += continuation_char;
        }
        result += body.substr(0, chunk);
        body = body.substr(chunk);
    }
    return result;
}

} // namespace fmig
