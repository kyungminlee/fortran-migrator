#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace fmig {

// A single source-level replacement: replace [offset, offset+length) with
// replacement_text.
struct Transformation {
    size_t offset;             // Byte offset in the source file.
    size_t length;             // Number of bytes to replace.
    std::string replacement;   // New text.
    std::string description;   // Human-readable label (for --dry-run).
};

// Validate a set of transformations: check for overlaps and sort by offset
// descending (so applying them in order doesn't shift later offsets).
// Returns false if overlapping transformations are found.
bool validate_and_sort(std::vector<Transformation>& transforms);

// Apply a sorted (descending offset) list of transformations to source text.
std::string apply_transformations(std::string_view source,
                                  const std::vector<Transformation>& transforms);

} // namespace fmig
