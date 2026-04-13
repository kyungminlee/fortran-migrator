#pragma once

#include <string>
#include <string_view>

namespace fmig {

// Detected source format of a Fortran file.
enum class FortranForm {
    Fixed, // .f, .for — columns matter
    Free,  // .f90, .f95, .F90
};

// Determine the source form from a file extension.
FortranForm detect_form(std::string_view extension);

// For fixed-form source, check whether a line is a comment.
bool is_fixed_form_comment(std::string_view line);

// For fixed-form source, check whether a line is a continuation.
bool is_fixed_form_continuation(std::string_view line);

// For fixed-form source, extract the statement portion (columns 7-72).
std::string_view fixed_form_statement(std::string_view line);

// Reformat a fixed-form line to stay within the 72-column limit,
// splitting into continuation lines if necessary.  Uses the given
// continuation character (default '+').
std::string reformat_fixed_line(std::string_view line,
                                char continuation_char = '+');

} // namespace fmig
