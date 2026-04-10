#pragma once

#include <filesystem>
#include <string>
#include <string_view>

#include "symbol_database.h"

namespace fmig {

// Compute the output filename for a migrated source file.
// e.g., "dgemm.f" → "qgemm.f" for KIND=16.
std::string target_filename(const std::string& original_name,
                            const SymbolDatabase& symbols,
                            int target_kind);

// Write the transformed source to an output file, creating parent
// directories as needed.
bool write_output(const std::filesystem::path& output_path,
                  std::string_view content);

} // namespace fmig
