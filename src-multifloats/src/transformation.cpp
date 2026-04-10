#include "transformation.h"

#include <algorithm>
#include <iostream>

namespace fmig {

bool validate_and_sort(std::vector<Transformation>& transforms) {
    // Sort by offset ascending first to check for overlaps.
    std::sort(transforms.begin(), transforms.end(),
              [](const Transformation& a, const Transformation& b) {
                  return a.offset < b.offset;
              });

    // Check for overlaps.
    for (size_t i = 1; i < transforms.size(); ++i) {
        if (transforms[i].offset < transforms[i - 1].offset + transforms[i - 1].length) {
            std::cerr << "Overlapping transformations at offsets "
                      << transforms[i - 1].offset << " and "
                      << transforms[i].offset << "\n";
            return false;
        }
    }

    // Re-sort descending so we can apply back-to-front without offset shifts.
    std::sort(transforms.begin(), transforms.end(),
              [](const Transformation& a, const Transformation& b) {
                  return a.offset > b.offset;
              });

    return true;
}

std::string apply_transformations(std::string_view source,
                                  const std::vector<Transformation>& transforms) {
    std::string result(source);
    // Transforms must be sorted descending by offset.
    for (const auto& t : transforms) {
        result.replace(t.offset, t.length, t.replacement);
    }
    return result;
}

} // namespace fmig
