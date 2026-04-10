#include "rewriter.h"

#include <algorithm>
#include <cctype>
#include <fstream>

namespace fmig {

std::string target_filename(const std::string& original_name,
                            const SymbolDatabase& symbols,
                            int target_kind) {
    // Extract the stem (without extension) and extension.
    std::filesystem::path p(original_name);
    std::string stem = p.stem().string();
    std::string ext = p.extension().string();

    // Try to look up the stem as a symbol to get its renamed form.
    std::string new_stem = symbols.target_name(stem, target_kind);
    if (new_stem.empty())
        return original_name; // No rename applicable.

    // Preserve the case of the original filename.
    // If original was lowercase, make the new stem lowercase too.
    bool all_lower = std::all_of(stem.begin(), stem.end(),
                                 [](unsigned char c) {
                                     return !std::isalpha(c) || std::islower(c);
                                 });
    if (all_lower) {
        std::transform(new_stem.begin(), new_stem.end(), new_stem.begin(),
                       [](unsigned char c) { return std::tolower(c); });
    }

    return new_stem + ext;
}

bool write_output(const std::filesystem::path& output_path,
                  std::string_view content) {
    std::filesystem::create_directories(output_path.parent_path());
    std::ofstream out(output_path, std::ios::binary);
    if (!out)
        return false;
    out.write(content.data(), static_cast<std::streamsize>(content.size()));
    return out.good();
}

} // namespace fmig
